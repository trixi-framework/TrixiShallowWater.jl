
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the SWE-Exner equations with a discontinuous sediment bed 
# to test well-balancedness

# Equations with Grass model
equations = ShallowWaterExnerEquations1D(gravity_constant = 9.81, H0 = 1.0,
                                         rho_f = 0.5, rho_s = 1.0, porosity = 0.5,
                                         friction = ManningFriction(n = 0.01),
                                         sediment_model = GrassModel(A_g = 0.01))

function initial_condition_steady_state(x, t,
                                        equations::ShallowWaterExnerEquations1D)
    hv = 0.0

    # discontinuous sediment bed
    if -2 < x[1] <= 0
        h = 0.25
        h_b = equations.threshold_limiter
    else
        h = equations.threshold_limiter
        h_b = 0.5
    end

    return SVector(h, hv, h_b)
end

initial_condition = initial_condition_steady_state

boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = water_sediment_height)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the TreeMesh

coordinates_min = -2.0
coordinates_max = 2.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solver
tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 5.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation
# use a Runge-Kutta method with CFL-based time step
sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);

summary_callback() # print the timer summary
