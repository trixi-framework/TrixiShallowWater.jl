
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

# Semidiscretization of the shallow water moment equations to test well-balancedness for a 
# lake-at-rest configuration with smooth bottom topography and two moments.
equations = ShallowWaterMomentEquations1D(gravity = 9.812, H0 = 1.75,
                                          n_moments = 2)

# For testing purposes this initial condition is also used together with the
# `ShallowWaterLinearizedMomentEquations1D`.
function initial_condition_well_balanced(x, t,
                                         equations::Union{ShallowWaterMomentEquations1D,
                                                          ShallowWaterLinearizedMomentEquations1D})
    # Initial lake-at-rest configuration
    H = 1.75
    v = 0.0
    a = zeros(equations.n_moments)

    # Set smooth bottom topography
    b = 1 / exp(0.5 * x[1]^2)

    h = H - b

    return SVector(h, h * v, (h * a)..., b)
end

initial_condition = initial_condition_well_balanced

###############################################################################
# Get the DG approximation space

volume_flux = (flux_careaga_etal, flux_nonconservative_careaga_etal)
surface_flux = (flux_careaga_etal, flux_nonconservative_careaga_etal)

solver = DGSEM(polydeg = 3,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Create the TreeMesh for the domain [-1, 1]

coordinates_min = -4.0
coordinates_max = 4.0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 5, # 2^refinement_level
                n_cells_max = 10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition_slip_wall)

###############################################################################
# ODE solver
tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks
summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (entropy,
                                                                 lake_at_rest_error))
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
