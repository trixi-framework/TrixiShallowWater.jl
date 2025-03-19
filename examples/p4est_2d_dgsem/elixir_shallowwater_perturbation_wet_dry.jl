
using Downloads: download
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

##
# we might need this change in the Trixi main
# @inline function velocity(u, equations::ShallowWaterEquations2D)
#     h, h_v1, h_v2, _ = u

# -    v1 = h_v1 / h
# -    v2 = h_v2 / h
# +    # v1 = h_v1 / h
# +    # v2 = h_v2 / h
# +    # Velocity desingularization
# +    v1 = (2.0 * h * h_v1) / (h^2 + max(h^2, 1e-6))
# +    v2 = (2.0 * h * h_v2) / (h^2 + max(h^2, 1e-6))
#     return SVector(v1, v2)
# end
#

# TODO: This elixir is for debugging purposes of the AMR + well-balanced mortars.
# Currently this crashes very early in the simulation. ONe possible explanation is that
# the projection / interpolation of the solution as it is refined / coarsened may need
# to use a similar strategy as the mortar projection where we ensure that a constant solution
# is projected on the dry elements. Further investigation is needed to see if this is indeed
# the case.

###############################################################################
# semidiscretization of the shallow water equations with a continuous
# bottom topography function

equations = ShallowWaterEquationsWetDry2D(gravity_constant = 9.812, H0 = 1.235)#, threshold_partially_wet=1e-4)#, threshold_limiter=1e-6) #, threshold_limiter=1e-4) #5*eps())

function initial_condition_perturbation(x, t, equations::ShallowWaterEquationsWetDry2D)
    # Calculate primitive variables
    H = equations.H0
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = (1.5 / exp(0.5 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
         +
         0.8 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
         -
         0.5 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))

    # for testing
    # b = zero(eltype(x))

    H = max(equations.H0, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_perturbation

# Wall BCs
boundary_condition = Dict(:all => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

# # Create the solver (for fully wet!)
# solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
#                volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

# Create the solver
basis = LobattoLegendreBasis(4)

# indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
#                                                      alpha_max = 1.0,
#                                                      alpha_min = 0.001,
#                                                      alpha_smooth = false,
#                                                      variable = Trixi.waterheight)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = Trixi.waterheight)
                                                    #  variable = waterheight_pressure)
                                                    #  variable = pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# This setup is for the curved, split form well-balancedness testing

# TODO: Make this mesh periodic to properly test for conservation
#       for this wet/dry testcase with AMR + well-balanced mortars

# Unstructured mesh with 24 cells of the square domain [-1, 1]^2
mesh_file = Trixi.download("https://gist.githubusercontent.com/efaulhaber/63ff2ea224409e55ee8423b3a33e316a/raw/7db58af7446d1479753ae718930741c47a3b79b7/square_unstructured_2.inp",
                           joinpath(@__DIR__, "square_unstructured_2.inp"))

# Affine type mapping to take the [-1,1]^2 domain from the mesh file
# and warp it as described in https://arxiv.org/abs/2012.12040
# Warping with the coefficient 0.2 is even more extreme.
function mapping_twist(xi, eta)
    y = eta + 0.15 * cos(1.5 * pi * xi) * cos(0.5 * pi * eta)
    x = xi + 0.15 * cos(0.5 * pi * xi) * cos(2.0 * pi * y)
    return SVector(x, y)
end

mesh = P4estMesh{2}(mesh_file, polydeg = 4,
                    mapping = mapping_twist,
                    initial_refinement_level = 0)

# Refine bottom left quadrant of each tree to level 2
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

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

function initial_condition_discontinuous_well_balancedness(x, t, element_id, equations::ShallowWaterEquationsWetDry2D)
    # Set the background values for velocity
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
    if element_id > 207 && element_id < 300 # for the forced hanging node grid with with < 3 in function above
        b = (0.75 / exp(0.5 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
            +
            0.4 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
            -
            0.25 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))
    end

    # Put in a discontinous purturbation using the element number
    if element_id in [232, 224, 225, 226, 227, 228, 229, 230]
        H = H + 1.3 # 1.8
    end

    # Clip the initialization to avoid negative water heights and division by zero
    h = max(equations.threshold_limiter, H - b)

    # Return the conservative variables
    return SVector(h, h * v1, h * v2, b)
end

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
  for j in eachnode(semi.solver), i in eachnode(semi.solver)
      x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations, semi.solver, i, j, element)
      u_node = initial_condition_discontinuous_well_balancedness(x_node, first(tspan), element, equations)
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

save_solution = SaveSolutionCallback(dt = 0.02, #interval=20,
                                     save_initial_solution = true,
                                     save_final_solution = true)

# # Define a better variable to use in the AMR indicator
# @inline function total_water_height(u, equations::ShallowWaterEquationsWetDry2D)
#   # return u[1]
#   return min(abs(u[1] + u[4] - equations.H0 - equations.threshold_limiter), abs(u[1] - equations.threshold_limiter)) + equations.H0
# end

# amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable = total_water_height),
#                                       base_level = 0,
#                                       med_level = 1, med_threshold = 1.1, # with adding back H0
#                                       max_level = 2, max_threshold = 1.3)

# amr_callback = AMRCallback(semi, amr_controller,
#                            interval = 1,
#                            adapt_initial_condition = false,
#                            adapt_initial_condition_only_refine = false)

stepsize_callback = StepsizeCallback(cfl = 0.8)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        # amr_callback,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(stage_limiter!),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            adaptive = false,
            save_everystep = false, callback = callbacks);
