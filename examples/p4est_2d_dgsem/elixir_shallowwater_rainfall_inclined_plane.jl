
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the shallow water equations with a rainfall and infiltration source term

equations = ShallowWaterEquations2D(gravity = 9.81)

# Initial condition with zero water height and a linearly inclined bottom topography to test
# rainfall and infiltration source terms. The setup for this test case is described in the paper:
#
# - J. Fernández-Pato, D. Caviedes-Voullième, P. García-Navarro (2016)
#   "Rainfall/runoff simulation with 2D full shallow water equations: Sensitivity analysis and
#   calibration of infiltration parameters"
#   [doi: 10.1016/j.jhydrol.2016.03.021](http://dx.doi.org/10.1016/j.jhydrol.2016.03.021)
function initial_condition_inclined_plane(x, t,
                                          equations::ShallowWaterEquations2D)
    b = 10.0 - 10.0 / 2000.0 * x[1]

    H = equations.threshold_limiter + b

    v1 = 0.0
    v2 = 0.0

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_inclined_plane

# Simple outflow boundary condition with fix for vanishing water heights.
function boundary_condition_outflow(u_inner, normal_direction::AbstractVector, x, t,
                                    surface_flux_functions,
                                    equations::ShallowWaterEquations2D)
    if u_inner[1] < 1e-8
        return boundary_condition_slip_wall(u_inner, normal_direction, x, t,
                                            surface_flux_functions, equations)
    else
        surface_flux_function, nonconservative_flux_function = surface_flux_functions
        # Impulse and bottom from inside, height from external state
        u_outer = SVector(1e-8, u_inner[2], u_inner[3], u_inner[4])

        # calculate the boundary flux
        flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction,
                                                     equations)

        return flux, noncons_flux
    end
end

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_wintermeyer_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(6)

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
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (0.0, 0.0)
coordinates_max = (2000.0, 20.0)

# Create P4estMesh
trees_per_dimension = (128, 2)
mesh = P4estMesh(trees_per_dimension, polydeg = 1,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 initial_refinement_level = 1,
                 periodicity = false)

boundary_conditions = (; x_neg = boundary_condition_slip_wall,
                       x_pos = boundary_condition_outflow,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

# Source term setup

# Prescribed time-dependent and spatially uniform rainfall pattern corresponding to Case 1-1 in 
# Section 3.2 of the reference:
# - J. Fernández-Pato, D. Caviedes-Voullième, P. García-Navarro (2016)
#  "Rainfall/runoff simulation with 2D full shallow water equations: Sensitivity analysis and
#  calibration of infiltration parameters"
#  [doi: 10.1016/j.jhydrol.2016.03.021](http://dx.doi.org/10.1016/j.jhydrol.2016.03.021)
function rain_pattern(x, t)
    if t <= 250 * 60
        return 1.25e-4
    else
        return 0.0
    end
end

# Manning friction source term
@inline function source_term_manning_friction(u, x, t,
                                              equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    n = 0.001 # friction coefficient
    h = (h^2 + max(h^2, 1e-8)) / (2 * h) # desingularization procedure

    ## Compute the common friction term
    Sf = -equations.gravity * n^2 * h^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    return SVector(zero(eltype(x)), Sf * hv_1, Sf * hv_2, zero(eltype(x)))
end

# Combined source term that includes rainfall / infiltration and bottom friction
function source_terms(u, x, t, equations::ShallowWaterEquations2D)
    infiltration_model = HortonModel(1.977e-4, 3.272e-5, 2.43e-3)
    precipitation_rate = rain_pattern
    src_rain = SourceTermsRain(precipitation_rate, infiltration_model, equations)(u, x, t,
                                                                                  equations)
    src_friction = source_term_manning_friction(u, x, t, equations)
    return src_rain + src_friction
end

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_conditions,
                                    source_terms = source_terms)
###############################################################################
# ODE solver

tspan = (0.0, 18000.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 10000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        save_solution)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(; stage_limiter! = stage_limiter!);
            dt = 0.1, ode_default_options()..., callback = callbacks, adaptive = false);
