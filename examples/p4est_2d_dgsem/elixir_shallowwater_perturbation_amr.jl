
using Downloads: download
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations with a continuous
# bottom topography function and a perturbation in the water height

equations = ShallowWaterEquationsWetDry2D(gravity_constant = 9.812, H0 = 2.1)

function initial_condition_perturbation(x, t, equations::ShallowWaterEquationsWetDry2D)
    # Calculate primitive variables
    H = equations.H0
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = (1.75 / exp(0.5 * ((x1 - 1.0)^2 + (x2 + 1.0)^2))
         +
         0.8 / exp(0.5 * ((x1 + 1.0)^2 + (x2 - 1.0)^2))
         -
         0.5 / exp(3.5 * ((x1 - 0.4)^2 + (x2 - 0.325)^2)))

    # Waterheight perturbation
    H = H + 0.5 * exp(-40.0 * ((x[1])^2 + (x[2])^2))

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_perturbation

boundary_condition = Dict(:x_neg => boundary_condition_slip_wall,
                          :y_neg => boundary_condition_slip_wall,
                          :x_pos => boundary_condition_slip_wall,
                          :y_pos => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

# Create the solver
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# This setup is for the curved, split form well-balancedness testing

# Affine type mapping to take the [-1,1]^2 domain from the mesh file
# and warp it as described in https://arxiv.org/abs/2012.12040
function mapping_twist(xi, eta)
    y = eta + 0.2 * cos(1.5 * pi * xi) * cos(0.5 * pi * eta)
    x = xi + 0.2 * cos(0.5 * pi * xi) * cos(2.0 * pi * y)
    return SVector(x, y)
end

# The mesh below can be made periodic if desired

# Create P4estMesh with 4 x 4 trees
trees_per_dimension = (4, 4)
mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 mapping = mapping_twist,
                 initial_refinement_level = 0,
                 periodicity = false)

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

tspan = (0.0, 3.0)
ode = semidiscretize(semi, tspan)

###############################################################################

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.04,
                                     save_initial_solution = true,
                                     save_final_solution = true)

# Define the perturbation of water height as a  variable to use in the AMR indicator
@inline function waterheight_total(u, equations::ShallowWaterEquationsWetDry2D)
    return u[1] + u[4]
end

amr_controller = ControllerThreeLevel(semi,
                                      IndicatorMax(semi, variable = waterheight_total),
                                      base_level = 1,
                                      med_level = 2, med_threshold = 2.01,
                                      max_level = 4, max_threshold = 2.15)

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 1,
                           adapt_initial_condition = false,
                           adapt_initial_condition_only_refine = false)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        amr_callback,
                        stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            adaptive = false,
            save_everystep = false, callback = callbacks);
