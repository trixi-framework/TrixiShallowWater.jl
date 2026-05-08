
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

# Semidiscretization of the shallow water moment equations with two moments
equations = ShallowWaterMomentEquations1D(gravity = 1.0, n_moments = 2)

# Initial condition with a smooth wave wave in a periodic domain.
# See section 4.2 in the paper:
#   Julian Koellermeier and Marvin Rominger (2020)
#   "Analysis and Numerical Simulation of Hyperbolic Shallow Water Moment Equations"
#   [DOI: 10.4208/cicp.OA-2019-0065](https://doi.org/10.4208/cicp.OA-2019-0065)
function initial_condition_smooth_periodic_wave(x, t,
                                                equations::Union{ShallowWaterMomentEquations1D,
                                                                 ShallowWaterLinearizedMomentEquations1D})
    H = 1.0 + exp(3.0 * cos(π * (x[1] + 0.5))) / exp(4.0)
    v = 0.25
    a = zeros(equations.n_moments)
    a[2] = -0.25

    b = 0.0

    return prim2cons(SVector(H, v, a..., b), equations)
end

initial_condition = initial_condition_smooth_periodic_wave

###############################################################################
# Get the DG approximation space

volume_flux = (flux_careaga_etal, flux_nonconservative_careaga_etal)
surface_flux = (FluxPlusDissipation(flux_careaga_etal,
                                    DissipationLaxFriedrichsEntropyVariables(Trixi.max_abs_speed)),
                flux_nonconservative_careaga_etal)

indicator_var(u, equations) = u[2]^3

basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = indicator_var)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux,)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Create the TreeMesh for the domain [-1, 1]

coordinates_min = -1.0
coordinates_max = 1.0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 6, # 2^refinement_level
                n_cells_max = 10_000,
                periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term_manning_friction;
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks
summary_callback = SummaryCallback()

analysis_interval = 200
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false)
alive_callback = AliveCallback(analysis_interval = analysis_interval)
stepsize_callback = StepsizeCallback(cfl = 0.9)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback)
###############################################################################
# run the simulation

sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,              # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()...,
            callback = callbacks,);
