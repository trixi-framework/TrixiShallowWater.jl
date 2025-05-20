
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

# Elixir useful for stress testing the AMR + wet/dry + well-balanced mortars.
# In the current state, AMR on wet/dry fronts remains sensitive to the limiting of the water height,
# desingularization of the velocities, and the design of the AMR indicator.
# Further, there remains an appreciable loss of conservation, with conservation errors
# of ≈1e-5 in the water height for this test case as the flow evolves.
# More investigation is needed to identify exactly why each aspect of this strange behavior occurs.
# For instance, the conservation loss  may be related to the choice of the CFL that guarantees
# a positive water height as well as a positive element-wise mean necessary in the limiting.

###############################################################################
# Semidiscretization of the multilayer shallow water equations with one layer
# to actually be the standard shallow water equations with a discontinuous
# bottom topography function for a perturbed water height
# on a nonconforming mesh with AMR

equations = ShallowWaterMultiLayerEquations2D(gravity = 9.812, H0 = 1.235,
                                              rhos = 1.0)

function initial_condition_perturbation(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Calculate primitive variables
    H = equations.H0
    v1 = zero(H)
    v2 = zero(H)

    x1, x2 = x
    b = (2.5 / exp(0.5 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
         +
         0.8 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
         -
         0.5 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations2D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    H = max(H, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_perturbation

# Wall BCs
boundary_condition = Dict(:all => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

# Create the solver
basis = LobattoLegendreBasis(4)

# Cannot simply use `waterheight` here for multilayer equations.
# Need a helper function to extract the relevant variable.
@inline function main_waterheight(u, equations)
    return waterheight(u, equations)[1]
end

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = main_waterheight)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the unstructured quad mesh from a file (downloads the file if not available locally)

# Unstructured mesh with 24 cells of the square domain [-1, 1]^2
mesh_file = Trixi.download("https://gist.githubusercontent.com/efaulhaber/63ff2ea224409e55ee8423b3a33e316a/raw/7db58af7446d1479753ae718930741c47a3b79b7/square_unstructured_2.inp",
                           joinpath(@__DIR__, "square_unstructured_2.inp"))

# Affine type mapping to take the [-1,1]^2 domain from the mesh file
# and warp it as described in https://arxiv.org/abs/2012.12040
# Warping with the coefficient 0.15 is even more extreme.
function mapping_twist(xi, eta)
    y = eta + 0.15 * cos(1.5 * pi * xi) * cos(0.5 * pi * eta)
    x = xi + 0.15 * cos(0.5 * pi * xi) * cos(2.0 * pi * y)
    return SVector(x, y)
end

mesh = P4estMesh{2}(mesh_file, polydeg = 4,
                    mapping = mapping_twist,
                    initial_refinement_level = 0)

# Refine bottom left quadrant of each tree to level 3
function refine_fn(p4est, which_tree, quadrant)
    quadrant_obj = unsafe_load(quadrant)
    if quadrant_obj.x == 0 && quadrant_obj.y == 0 && quadrant_obj.level < 3
        # return true (refine)
        return Cint(1)
    else
        # return false (don't refine)
        return Cint(0)
    end
end

# Refine recursively until each bottom left quadrant of a tree has level 3
# The mesh will be rebalanced before the simulation starts
refine_fn_c = @cfunction(refine_fn, Cint,
                         (Ptr{Trixi.p4est_t}, Ptr{Trixi.p4est_topidx_t},
                          Ptr{Trixi.p4est_quadrant_t}))
Trixi.refine_p4est!(mesh.p4est, true, refine_fn_c, C_NULL)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks, etc.

tspan = (0.0, 5.0)
ode = semidiscretize(semi, tspan)

function initial_condition_discontinuous_perturbation(x, t, element_id,
                                                      equations::ShallowWaterMultiLayerEquations2D)
    # Set the background values for velocity
    H = equations.H0
    v1 = zero(H)
    v2 = zero(H)

    x1, x2 = x

    # For the mesh file version of the testing
    b = (1.75 / exp(0.6 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
         +
         0.8 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
         -
         0.5 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))

    # Setup a discontinuous bottom topography using the element id number
    # Note that this requires to have the mesh exactly as created above,
    # any additional refinement changes the initial condition,
    IDs = [collect(114:133); collect(138:141); collect(156:164); collect(208:300)]
    if element_id in IDs
        b = (0.75 / exp(0.5 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
             +
             0.4 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
             -
             0.25 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))
    end

    # Put in a discontinuous perturbation using the element number
    if element_id in [232, 224, 225, 226, 227, 228, 229, 230]
        H = H + 1.6
    end

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations2D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    H = max(H, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for j in eachnode(semi.solver), i in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, j, element)
        u_node = initial_condition_discontinuous_perturbation(x_node, first(tspan), element,
                                                              equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
    end
end

###############################################################################

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.2,
                                     save_initial_solution = true,
                                     save_final_solution = true)

amr_controller = ControllerThreeLevel(semi,
                                      IndicatorMax(semi, variable = main_waterheight),
                                      base_level = 0,
                                      med_level = 1, med_threshold = 0.3,
                                      max_level = 4, max_threshold = 1.15)

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 5,
                           adapt_initial_condition = false,
                           adapt_initial_condition_only_refine = false)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        amr_callback,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(stage_limiter!),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            adaptive = false,
            save_everystep = false, callback = callbacks);
