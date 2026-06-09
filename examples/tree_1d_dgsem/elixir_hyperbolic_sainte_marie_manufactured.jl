using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the hyperbolic Sainte-Marie equations to test convergence
# against a smooth, periodic manufactured solution.

equations = HyperbolicSainteMarieEquations1D(gravity = 1.0, h_ref = 2.0, alpha = 3.0)

initial_condition = initial_condition_manufactured

###############################################################################
# Get the DG approximation space
alpha_coefficients = (1 / 2, 1.0, 2 / 3)
volume_flux = (FluxArtianoEtal(alpha_coefficients...),
               FluxNonConservativeArtianoEtal(alpha_coefficients...))
surface_flux = (FluxPlusDissipation(FluxArtianoEtal(alpha_coefficients...),
                                    DissipationLocalLaxFriedrichs(Trixi.max_abs_speed)),
                FluxNonConservativeArtianoEtal(alpha_coefficients...))

solver = DGSEM(polydeg = 3,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = 0.0
coordinates_max = 1.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semidiscretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_manufactured,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver

tspan = (0.0, 0.1)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation

abstol = 1.0e-9
reltol = 1.0e-9
sol = solve(ode, SSPRK43(thread = Trixi.Threaded());
            abstol = abstol, reltol = reltol,
            save_everystep = false, callback = callbacks);
