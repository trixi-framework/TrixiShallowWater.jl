
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the shallow water equations

# By passing only a single value for rhos, the system recovers the standard shallow water equations
equations = ShallowWaterMultiLayerEquations1D(gravity = 9.81, rhos = 1.0)

"""
    initial_condition_parabolic_bowl(x, t, equations:: ShallowWaterMultiLayerEquations1D)

Well-known initial condition to test the [`hydrostatic_reconstruction_ersing_etal`](@ref) and its
wet-dry mechanics. This test has analytical solutions. The initial condition is defined by the
analytical solution at time t=0. The bottom topography defines a bowl and the water level is given
by an oscillating lake.

The original test and its analytical solution in two dimensions were first presented in
- William C. Thacker (1981)
  Some exact solutions to the nonlinear shallow-water wave equations
  [DOI: 10.1017/S0022112081001882](https://doi.org/10.1017/S0022112081001882).

The particular setup below is taken from Section 6.2 of
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and Timothy Warburton (2018)
  An entropy stable discontinuous Galerkin method for the shallow water equations on
  curvilinear meshes with wet/dry fronts accelerated by GPUs
  [DOI: 10.1016/j.jcp.2018.08.038](https://doi.org/10.1016/j.jcp.2018.08.038).
"""
function initial_condition_parabolic_bowl(x, t,
                                          equations::ShallowWaterMultiLayerEquations1D)
    a = 1
    h_0 = 0.1
    sigma = 0.5
    ω = sqrt(2 * equations.gravity * h_0) / a

    v = -sigma * ω * sin(ω * t)

    b = h_0 * x[1]^2 / a^2

    H = sigma * h_0 / a^2 * (2 * x[1] * cos(ω * t) - sigma) + h_0

    #=
    It is mandatory to shift the water level at dry areas to make sure the water height h
    stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    for the ShallowWaterMultiLayerEquations1D and added to the initial condition if h = 0.
    This default value can be changed within the constructor call depending on the simulation setup.
    =#
    H = max(H, b + equations.threshold_limiter)
    return prim2cons(SVector(H, v, b), equations)
end

initial_condition = initial_condition_parabolic_bowl

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
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs(max_abs_speed_naive)),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(5)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Create the TreeMesh for the domain [-2, 2]

coordinates_min = -2.0
coordinates_max = 2.0

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 6,
                n_cells_max = 10_000)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (energy_kinetic,
                                                                 energy_internal))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks);
