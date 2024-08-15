
using Downloads: download
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations with a continuous
# bottom topography function

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

# Wall BCs
boundary_condition = Dict(:all => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

# Create the solver
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# This setup is for the curved, split form well-balancedness testing

# Unstructured mesh with 24 cells of the square domain [-1, 1]^2
mesh_file = Trixi.download("https://gist.githubusercontent.com/efaulhaber/63ff2ea224409e55ee8423b3a33e316a/raw/7db58af7446d1479753ae718930741c47a3b79b7/square_unstructured_2.inp",
                           joinpath(@__DIR__, "square_unstructured_2.inp"))

# Affine type mapping to take the [-1,1]^2 domain from the mesh file
# and warp it as described in https://arxiv.org/abs/2012.12040
# Warping with the coefficient 0.2 is even more extreme.
function mapping_twist(xi, eta)
    y = eta + 0.175 * cos(1.5 * pi * xi) * cos(0.5 * pi * eta)
    x = xi + 0.175 * cos(0.5 * pi * xi) * cos(2.0 * pi * y)
    return SVector(x, y)
end

mesh = P4estMesh{2}(mesh_file, polydeg = 3,
                    mapping = mapping_twist,
                    initial_refinement_level = 0)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks, etc.

tspan = (0.0, 0.5)
ode = semidiscretize(semi, tspan)

###############################################################################

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.04,
                                     save_initial_solution = true,
                                     save_final_solution = true)

# Define a better variable to use in the AMR indicator
@inline function total_water_height(u, equations::ShallowWaterEquationsWetDry2D)
  return u[1] + u[4]
end

amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable = total_water_height),
                                      base_level = 0,
                                      med_level = 1, med_threshold = 2.01,
                                      max_level = 4, max_threshold = 2.15)

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 1,
                           adapt_initial_condition = false, # to setup discontinuous bottom easier
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
            # adaptive = false,
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary