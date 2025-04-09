
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the SWE-Exner equations

# Equations with Grass model
equations = ShallowWaterExnerEquations1D(gravity_constant = 9.812, H0 = 1.0,
                                         rho_f = 0.5, rho_s = 1.0, porosity = 0.5,
                                         friction = ManningFriction(n = 0.0),
                                         sediment_model = GrassModel(A_g = 0.001),
                                         threshold_desingularization = 1e-6)

"""
    initial_condition_beach(x, t, equations:: ShallowWaterExnerEquations1D)

Initial condition to simulate a wave running towards a beach and crashing. Difficult test
including both wetting and drying in the domain using slip wall boundary conditions.
The bottom topography is altered to be differentiable on the domain [0,8] and
differs from the reference below.

The water height and speed functions used here, are adapted from the initial condition
found in section 5.2 of the paper:
  - Andreas Bollermann, Sebastian Noelle, Maria Lukáčová-Medvid’ová (2011)
    Finite volume evolution Galerkin methods for the shallow water equations with dry beds\n
    [DOI: 10.4208/cicp.220210.020710a](https://dx.doi.org/10.4208/cicp.220210.020710a)
"""
function initial_condition_beach(x, t, equations::ShallowWaterExnerEquations1D)
    D = 1
    delta = 0.02
    gamma = sqrt((3 * delta) / (4 * D))
    x_a = sqrt((4 * D) / (3 * delta)) * acosh(sqrt(20))

    f = D + 40 * delta * sech(gamma * (8 * x[1] - x_a))^2

    # steep curved beach
    h_b = 0.01 + 99 / 409600 * 4^x[1]

    if x[1] >= 6
        H = h_b
        v = 0.0
    else
        H = f
        v = sqrt(equations.gravity / D) * H
    end

    H = max(H, h_b + equations.threshold_limiter)
    return prim2cons(SVector(H, v, h_b), equations)
end

initial_condition = initial_condition_beach
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxPlusDissipation(flux_ersing_etal, DissipationLocalLaxFriedrichs()),
                flux_nonconservative_ersing_etal,
                hydrostatic_reconstruction_ersing_etal)

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
# Create the TreeMesh for the domain [0, 8]

coordinates_min = 0.0
coordinates_max = 8.0

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.2,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation
# use a Runge-Kutta method with CFL-based time step
sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);

summary_callback() # print the timer summary
