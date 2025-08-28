using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations with three layers

equations = ShallowWaterMultiLayerEquations2D(gravity = 1.1,
                                              rhos = (0.9, 1.0, 1.1))

initial_condition = initial_condition_convergence_test

###############################################################################
# Get the DG approximation space

polydeg = 3
volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)

basis = LobattoLegendreBasis(polydeg)
limiter_idp = SubcellLimiterIDP(equations, basis;)
volume_integral = VolumeIntegralSubcellLimiting(limiter_idp;
                                                volume_flux_dg = volume_flux,
                                                volume_flux_fv = surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the P4estMesh and setup a periodic mesh

coordinates_min = (0.0, 0.0) # minimum coordinates (min(x), min(y))
coordinates_max = (sqrt(2.0), sqrt(2.0))  # maximum coordinates (max(x), max(y))

# Create P4estMesh with 4 x 4 elements
trees_per_dimension = (2, 2)
mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 initial_refinement_level = 1,
                 periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 500,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.25)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

stage_callbacks = (SubcellLimiterIDPCorrection(),)

sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = stage_callbacks);
                  dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                  ode_default_options()...,
                  callback = callbacks, adaptive = false);
