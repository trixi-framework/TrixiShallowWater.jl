
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations to test inflow/outflow boundary conditions

equations = ShallowWaterEquations2D(gravity = 9.81)

# Setup initial conditions for a smooth channel flow with constant water height and velocity
function initial_condition_channel_flow(x, t, equations::ShallowWaterEquations2D)
    H = 1.0
    v1 = 0.4
    v2 = 0.4
    b = 0.0

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_channel_flow

# Setup boundary conditions.
# At the inlet, we prescribe constant momentum to simulate inflow. 
# At the outlet, we prescribe the water height as a time-dependent cosine wave.
boundary_condition_inflow = BoundaryConditionMomentum(0.4, 0.4, equations)
boundary_condition_outflow = BoundaryConditionWaterHeight(t -> 1.0 + 0.1 * cos(π / 2 * t),
                                                          equations)

boundary_conditions = Dict(:Bottom => boundary_condition_inflow,
                           :Top => boundary_condition_outflow,
                           :Left => boundary_condition_slip_wall,
                           :Right => boundary_condition_slip_wall)

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
# Get the unstructured quad mesh from a file (downloads the file if not available locally)
default_mesh_file = joinpath(@__DIR__, "mesh_wobbly_channel.mesh")
isfile(default_mesh_file) ||
    Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/431baa423ce86aadba70d38a3194947b/raw/50914cb30e72e9a58d4723e161476435c6dea182/mesh_wobbly_channel.mesh",
                   default_mesh_file)
mesh_file = default_mesh_file

mesh = UnstructuredMesh2D(mesh_file, periodicity = false)

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

save_solution = SaveSolutionCallback(dt = 0.4,
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
