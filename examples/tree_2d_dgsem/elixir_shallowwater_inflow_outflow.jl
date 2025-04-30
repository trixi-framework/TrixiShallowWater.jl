using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations to test inflow/outflow boundary conditions

equations = ShallowWaterEquationsWetDry2D(gravity = 9.81)

# Setup initial conditions for a smooth channel flow with constant water height and velocity
function initial_condition_channel_flow(x, t, equations::ShallowWaterEquationsWetDry2D)
    H = 1.0
    v1 = -0.1
    v2 = -0.1
    b = 0.0

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_channel_flow

# Setup boundary conditions.
# At the inlet, we prescribe the momentum, starting with a negative value to simulate outflow,
# and gradually transitioning to a positive value to simulate inflow. At the outlet, we prescribe
# the water height as a time-dependent cosine wave. This setup is designed to test the behavior
# of the boundary conditions under both inflow and outflow scenarios.
boundary_condition_inflow = BoundaryConditionMomentum(t -> -0.1 + 0.05 * t,
                                                      t -> -0.1 + 0.05 * t,
                                                      equations)
boundary_condition_outflow = BoundaryConditionWaterHeight(t -> 1.0 + 0.1 * cos(π / 2 * t),
                                                          equations)
boundary_conditions = (x_neg = boundary_condition_inflow,
                       x_pos = boundary_condition_outflow,
                       y_neg = boundary_condition_inflow,
                       y_pos = boundary_condition_outflow)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxPlusDissipation(flux_wintermeyer_etal, DissipationLocalLaxFriedrichs()),
                flux_nonconservative_wintermeyer_etal)
polydeg = 3
solver = DGSEM(polydeg = polydeg,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Create a non-periodic TreeMesh
coordinates_min = (-10.0, -10.0)
coordinates_max = (10.0, 10.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000, periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solver

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 1.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
