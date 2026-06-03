
OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water Exner equations with Meyer-Peter-Mueller
# sediment closure and discontinuous initial data

# Academic test case of entropy conservation.
# The errors from the analysis callback are not important but `∑∂S/∂U ⋅ Uₜ` is.
# If the Manning coefficient `n = 0`, then `∑∂S/∂U ⋅ Uₜ` should be around machine roundoff.
# If the Manning coefficient `n > 0`, then `∑∂S/∂U ⋅ Uₜ` should be negative.

# Equations with Meyer-Peter-Mueller model
equations = ShallowWaterExnerEquations2D(gravity = 10.0, rho_f = 0.5,
                                         rho_s = 1.0, porosity = 0.5,
                                         friction = ManningFriction(n = 0.01),
                                         sediment_model = MeyerPeterMueller(theta_c = 0.0,
                                                                            d_s = 1e-3))

# Note, this initial condition is used to compute errors in the analysis callback but the initialization is
# overwritten by `initial_condition_ec_discontinuous_bottom` below.
function Trixi.initial_condition_weak_blast_wave(x, t,
                                                 equations::ShallowWaterExnerEquations2D)
    # Set up polar coordinates
    RealT = eltype(x)
    inicenter = SVector(convert(RealT, 0.7), convert(RealT, 0.7))
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    phi = atan(y_norm, x_norm)
    sin_phi, cos_phi = sincos(phi)

    # Calculate primitive variables
    H = r > 0.5f0 ? 3.25f0 : 4.0f0
    v1 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * cos_phi
    v2 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * sin_phi
    h_b = 0 # by default assume there is no bottom topography

    return prim2cons(SVector(H, v1, v2, h_b), equations)
end

initial_condition = initial_condition_weak_blast_wave

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

solver = DGSEM(polydeg = 4,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 2,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term_bottom_friction,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Workaround to set a discontinuous bottom topography and initial condition for debugging and testing.

# Alternative version of the initial condition used to setup a truly discontinuous
# bottom topography function and initial condition.
# In contrast to the usual signature of initial conditions, this one get passed the
# `element_id` explicitly. In particular, this initial conditions works as intended
# only for the TreeMesh2D with initial_refinement_level=2.
function initial_condition_ec_discontinuous_bottom(x, t, element_id,
                                                   equations::ShallowWaterExnerEquations2D)
    # Set up polar coordinates
    inicenter = SVector(0.7, 0.7)
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    phi = atan(y_norm, x_norm)
    sin_phi, cos_phi = sincos(phi)

    # Set the background values
    H = 4.25
    v1 = 0.0
    v2 = 0.0
    h_b = 0.0

    # setup the discontinuous water height and velocities
    if element_id == 10
        H = 5.0
        v1 = 0.1882 * cos_phi
        v2 = 0.1882 * sin_phi
    end

    # Setup a discontinuous bottom topography using the element id number
    if element_id == 7
        h_b = 2.0 + 0.5 * sin(2.0 * pi * x[1]) + 0.5 * cos(2.0 * pi * x[2])
    end

    return prim2cons(SVector(H, v1, v2, h_b), equations)
end

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for j in eachnode(semi.solver), i in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, j, element)
        u_node = initial_condition_ec_discontinuous_bottom(x_node, first(tspan), element,
                                                           equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
    end
end

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.5,
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
