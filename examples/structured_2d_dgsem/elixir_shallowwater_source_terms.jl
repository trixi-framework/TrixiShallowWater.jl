
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations

equations = ShallowWaterEquations2D(gravity = 9.81)

initial_condition = initial_condition_convergence_test

# Up to Trixi.jl version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `const flux_lax_friedrichs = FluxLaxFriedrichs(), i.e., `FluxLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 of Trixi.jl the default wave speed used in the LLF flux to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `FluxLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.
surface_flux = (FluxLaxFriedrichs(max_abs_speed_naive), flux_nonconservative_fjordholm_etal)
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg = 3,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (0.0, 0.0)
coordinates_max = (sqrt(2.0), sqrt(2.0))

cells_per_dimension = (8, 8)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
