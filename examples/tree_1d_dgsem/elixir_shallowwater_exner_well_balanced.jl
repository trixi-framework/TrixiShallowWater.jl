
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the SWE-Exner equations with a discontinuous sediment bed
# to test well-balancedness

# Equations with Meyer-Peter-Mueller sedimentation model
equations = ShallowWaterExnerEquations1D(gravity = 10.0, H0 = 1.0,
                                         rho_f = 0.5, rho_s = 1.0, porosity = 0.5,
                                         friction = ManningFriction(n = 0.01),
                                         sediment_model = MeyerPeterMueller(theta_c = 0.047,
                                                                            d_s = 1e-3))

function initial_condition_steady_state(x, t,
                                        equations::ShallowWaterExnerEquations1D)
    hv = 0.0

    # discontinuous sediment bed
    if -2.0 < x[1] < 0.0
        h_s = 0.3
    else
        h_s = 0.1
    end

    h = 1.0 - h_s

    return SVector(h, hv, h_s)
end

initial_condition = initial_condition_steady_state

boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
# Up to Trixi.jl version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `const flux_lax_friedrichs = FluxLaxFriedrichs(), i.e., `FluxLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 the default wave speed used in the LLF flux to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `FluxLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.
surface_flux = (FluxPlusDissipation(flux_ersing_etal,
                                    DissipationLocalLaxFriedrichs(max_abs_speed_naive)),
                flux_nonconservative_ersing_etal)

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassner(equations, basis,
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
                                    source_terms = source_term_bottom_friction,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solver

tspan = (0.0, 200.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
