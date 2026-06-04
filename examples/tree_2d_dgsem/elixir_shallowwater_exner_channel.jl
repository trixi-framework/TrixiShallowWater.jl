
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the shallow water Exner equations for a channel flow problem
# with sediment transport

equations = ShallowWaterExnerEquations2D(gravity = 9.81, H0 = 10.0,
                                         rho_f = 1.0, rho_s = 1.0,
                                         porosity = 0.4,
                                         sediment_model = GrassModel(A_g = 0.01))

# Initial condition for a channel flow problem over a sediment hump.
# Discussed at length in the first reference below including how to estimate
# the spreading angle of the dune as it evolves into a star shape.
# The second reference contains numerical results for comparison in section 4.3.
# - H. J. de Vriend (1987)
#   2DH mathematical modelling of morphological evolutions in shallow water
#   [DOI: 10.1016/0378-3839(87)90037-8](https://doi.org/10.1016/0378-3839(87)90037-8)
# - F. Benkhaldoun, S. Sahmim, M. Seaïd (2010)
#   A two-dimensional finite volume morphodynamic model on unstructured triangular grids
#   [DOI: 10.1002/fld.2129](https://doi.org/10.1002/fld.2129)
function initial_condition_channel(x, t,
                                   equations::ShallowWaterExnerEquations2D)
    # Compute the sediment and water height
    h_b = zero(eltype(x))
    if 300 <= x[1] <= 500 && 400 <= x[2] <= 600
        h_b += (sin(pi * (x[1] - 300) / 200))^2 * (sin(pi * (x[2] - 400) / 200))^2
    end
    h = equations.H0 - h_b

    # Set the background values for momenta
    hv1 = 10
    hv2 = zero(eltype(x))

    return SVector(h, hv1, hv2, h_b)
end

initial_condition = initial_condition_channel

function boundary_condition_subcritical_outflow(u_inner, orientation, direction,
                                                x, t,
                                                surface_flux_functions,
                                                equations::ShallowWaterExnerEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # Impulse from inside, height and bottom from outside
    # This acts as a weakly nonreflective BC, but can be improved
    u_boundary = SVector(equations.H0, u_inner[2], u_inner[3], zero(eltype(u_inner)))

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
                                                     equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
                                                     equations)
    end

    return flux, noncons_flux
end

boundary_conditions = (; x_neg = BoundaryConditionDirichlet(initial_condition_channel),
                       x_pos = boundary_condition_subcritical_outflow,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxPlusDissipation(flux_ersing_etal, dissipation_roe),
                flux_nonconservative_ersing_etal)
solver = DGSEM(polydeg = 4,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh

coordinates_min = (0.0, 0.0)
coordinates_max = (1000.0, 1000.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solver

tspan = (0.0, 36_000.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

save_solution = SaveSolutionCallback(dt = 300.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback, save_solution)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
