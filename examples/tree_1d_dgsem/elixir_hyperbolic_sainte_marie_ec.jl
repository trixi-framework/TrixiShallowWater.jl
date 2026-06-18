using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the hyperbolic Sainte-Marie equations
# for a smooth and periodic initial condition to test entropy conservation

equations = HyperbolicSainteMarieEquations1D(gravity = 1.0, h_ref = 0.1)

function initial_condition_periodic(x, t, equations::HyperbolicSainteMarieEquations1D)
    h = 1 + exp(sinpi(2 * x[1]))
    v = 1
    w = 1
    p = 10
    b = 0.1 * exp(sinpi(2 * x[1]))
    H = h + b
    return prim2cons(SVector(H, v, w, p, b), equations)
end

initial_condition = initial_condition_periodic
###############################################################################
# Get the DG approximation space
alpha_coefficients = (1 / 2, 1.0, 2 / 3)
volume_flux = (FluxArtianoEtal(alpha_coefficients...),
               FluxNonConservativeArtianoEtal(alpha_coefficients...))
solver = DGSEM(polydeg = 3, surface_flux = volume_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
boundary_condition = BoundaryConditionDirichlet(initial_condition)
###############################################################################
# Get the TreeMesh and setup a periodic mesh
coordinates_min = 0.0
coordinates_max = 1.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term_hyperbolic_sainte_marie,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver
tspan = (0.0, 0.2)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (entropy,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback)

###############################################################################
# run the simulation

abstol = 1.0e-6
reltol = 1.0e-6
sol = solve(ode, SSPRK43(thread = Trixi.Threaded());
            abstol = abstol, reltol = reltol,
            save_everystep = false, callback = callbacks);
