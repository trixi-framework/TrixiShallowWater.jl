
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

# Stress test case for well-balancedness with wet/dry transition elements
# on a nonconforming mesh with a discontinuous bottom topography where
# this elixir allows for experimentation for well-balancedness with this complex test case.

# In the current state, recovering a well-balancedness error up to double precision machine
# round-off error relies on a special projection procedure at wet/dry fronts.
# However, there are no conservation issues and these conservation errors remain around
# double precision unit round-off for all time.

###############################################################################
# Semidiscretization of the multilayer shallow water equations with one layer
# to fallback to the standard shallow water equations.

equations = ShallowWaterMultiLayerEquations2D(gravity = 9.812, H0 = 1.235,
                                              rhos = 1.0)

# An initial condition with constant total water height and zero velocities to test well-balancedness.
# Note, this routine is used to compute errors in the analysis callback but the initialization is
# overwritten by `initial_condition_discontinuous_well_balancedness` below.
function initial_condition_well_balancedness(x, t,
                                             equations::ShallowWaterMultiLayerEquations2D)

    # Set the background values
    H = equations.H0
    v1 = zero(H)
    v2 = zero(H)

    x1, x2 = x
    b = (1.5 / exp(0.5 * ((x1 - 1)^2 + (x2 + 1)^2))
         +
         0.8 / exp(0.5 * ((x1 + 1)^2 + (x2 - 1)^2))
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

initial_condition = initial_condition_well_balancedness

boundary_condition = Dict(:all => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

# Up to Trixi.jl version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `DissipationLocalLaxFriedrichs(), i.e., `DissipationLocalLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 of Trixi.jl the default wave speed used in the LLF flux and dissipation operator to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `DissipationLocalLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs(max_abs_speed_naive)),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

# Create the solver
basis = LobattoLegendreBasis(3)

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

# Unstructured mesh with 24 cells of the square domain [-1, 1]^n
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

mesh = P4estMesh{2}(mesh_file, polydeg = 3,
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

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Workaround to set a discontinuous bottom topography for debugging and testing.
# The errors from the analysis callback are not important but the error for this lake at rest test case
# `∑|H0-(h+b)|` should be around machine roundoff
# In contrast to the usual signature of initial conditions, this one get passed the
# `element_id` explicitly. In particular, this initial conditions works as intended
# only for the specific mesh loaded above!
function initial_condition_discontinuous_well_balancedness(x, t, element_id,
                                                           equations::ShallowWaterMultiLayerEquations2D)
    # Set the background values
    H = equations.H0
    v1 = zero(H)
    v2 = zero(H)

    x1, x2 = x
    b = (1.75 / exp(0.6 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
         +
         0.8 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
         -
         0.5 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))

    # Setup a discontinuous bottom topography using the element id number
    IDs = [collect(114:133); collect(138:141); collect(156:164); collect(208:300)]
    if element_id in IDs
        b = (0.75 / exp(0.5 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
             +
             0.4 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
             -
             0.25 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))
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
        u_node = initial_condition_discontinuous_well_balancedness(x_node, first(tspan),
                                                                   element, equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
    end
end
###############################################################################

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 10.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(stage_limiter!),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            adaptive = false,
            save_everystep = false, callback = callbacks);
