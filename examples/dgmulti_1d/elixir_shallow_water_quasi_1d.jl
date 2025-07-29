using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the quasi 1d shallow water equations
# See Chan et al.  https://doi.org/10.48550/arXiv.2307.12089 for details

equations = ShallowWaterEquationsQuasi1D(gravity = 9.81)

initial_condition = initial_condition_convergence_test

volume_flux = (flux_chan_etal, flux_nonconservative_chan_etal)
# Up to Trixi.jl version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `DissipationLocalLaxFriedrichs(), i.e., `DissipationLocalLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 of Trixi.jl the default wave speed used in the LLF flux and dissipation operator to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `DissipationLocalLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.
surface_flux = (FluxPlusDissipation(flux_chan_etal,
                                    DissipationLocalLaxFriedrichs(max_abs_speed_naive)),
                flux_nonconservative_chan_etal)

dg = DGMulti(polydeg = 4, element_type = Line(), approximation_type = SBP(),
             surface_integral = SurfaceIntegralWeakForm(surface_flux),
             volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

cells_per_dimension = (8,)
mesh = DGMultiMesh(dg, cells_per_dimension,
                   coordinates_min = (0.0,), coordinates_max = (sqrt(2),),
                   periodicity = true)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg;
                                    source_terms = source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval, uEltype = real(dg))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation

sol = solve(ode, RDPK3SpFSAL49(); abstol = 1.0e-8, reltol = 1.0e-8,
            ode_default_options()..., callback = callbacks)
