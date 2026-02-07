
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations with a single layer and a bottom
# topography function for a blast wave test with a wet/dry front with discontinuous initial conditions.
# This example combines the subcell limiter with a velocity desingularization callback to ensure
# robustness at wet/dry fronts.

equations = ShallowWaterMultiLayerEquations2D(gravity = 9.81, H0 = 0.45,
                                              rhos = (1.0),
                                              threshold_desingularization = 1e-6)

function initial_condition_blast_wave(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Set up polar coordinates
    RealT = eltype(x)
    inicenter = SVector(convert(RealT, 2.0), convert(RealT, 2.0))
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    phi = atan(y_norm, x_norm)
    sin_phi, cos_phi = sincos(phi)

    # Calculate primitive variables
    H = r > 0.5f0 ? 2.0f0 : 4.0f0
    v1 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * cos_phi
    v2 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * sin_phi

    b = min(sqrt((2 - x[1])^2 + (2 - x[2])^2)^2, 4.0)

    H = max(H, b + equations.threshold_limiter)
    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_blast_wave

boundary_condition = (; x_neg = boundary_condition_slip_wall,
                      y_neg = boundary_condition_slip_wall,
                      x_pos = boundary_condition_slip_wall,
                      y_pos = boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

polydeg = 4
volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal_local_jump,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(polydeg)
limiter_idp = SubcellLimiterIDP(equations, basis;
                                positivity_variables_cons = ["h1"],
                                local_twosided_variables_cons = ["h1"],
                                positivity_correction_factor = 0.5,)
volume_integral = VolumeIntegralSubcellLimiting(limiter_idp;
                                                volume_flux_dg = volume_flux,
                                                volume_flux_fv = surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the P4estMesh

coordinates_min = (0.0, 0.0) # minimum coordinates (min(x), min(y))
coordinates_max = (4.0, 4.0) # maximum coordinates (max(x), max(y))

trees_per_dimension = (1, 1)

mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 initial_refinement_level = 5,
                 periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solver

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_polydeg = polydeg,
                                     extra_analysis_errors = (:conservation_error,))

stepsize_callback = StepsizeCallback(cfl = 0.25)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     extra_node_variables = (:limiting_coefficient,))

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################

# Setup the stage callbacks
stage_limiter! = VelocityDesingularization()
stage_callbacks = (SubcellLimiterIDPCorrection(), stage_limiter!)

# run the simulation
sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = stage_callbacks);
                  dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                  ode_default_options()...,
                  callback = callbacks, adaptive = false);
