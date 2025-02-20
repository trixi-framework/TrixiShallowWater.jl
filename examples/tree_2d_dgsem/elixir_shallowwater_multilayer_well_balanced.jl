
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations with a bottom topography function 
# to test well-balancedness

equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81, H0 = 0.6,
                                              rhos = (0.7, 0.8, 0.9, 1.0))

# An initial condition with constant total water height, zero velocities and a bottom topography to 
# test well-balancedness
function initial_condition_well_balanced(x, t, equations::ShallowWaterMultiLayerEquations2D)
    H = SVector(0.6, 0.55, 0.5, 0.45)
    v1 = zero(H)
    v2 = zero(H)
    b = (((x[1] - 0.5)^2 + (x[2] - 0.5)^2) < 0.04 ?
         0.2 * (cos(4 * pi * sqrt((x[1] - 0.5)^2 + (x[2] -
                                                    0.5)^2)) + 1) : 0.0)

    return prim2cons(SVector(H..., v1..., v2..., b),
                     equations)
end

initial_condition = initial_condition_well_balanced

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (0.0, 0.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solver

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

stepsize_callback = StepsizeCallback(cfl = 1.0)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
summary_callback() # print the timer summary
