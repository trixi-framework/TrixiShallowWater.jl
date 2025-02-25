
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break test with a
# discontinuous bottom topography function to test entropy conservation

equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 1.0,
                                              rhos = (0.9, 0.95, 1.0))

# This academic testcase sets up a discontinuous bottom topography
# function and initial condition to test entropy conservation.

function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 0.3 * exp(-0.5 * ((x[1])^2 + (x[2])^2))

    if x[1] < 0.0
        H = SVector(1.0, 0.8, 0.6)
    else
        H = SVector(0.9, 0.7, 0.5)
        b += 0.1
    end

    v1 = zero(H)
    v2 = zero(H)
    return prim2cons(SVector(H..., v1..., v2..., b),
                     equations)
end

initial_condition = initial_condition_dam_break

boundary_conditions = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a non-periodic mesh with wall boundary conditions

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)
###############################################################################
# ODE solver

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)
###############################################################################
#=
Workaround for TreeMesh2D to set true discontinuities for debugging and testing.
Essentially, this is a slight augmentation of the `compute_coefficients` where the `x` node values
passed here are slightly perturbed in order to set a true discontinuity that avoids the doubled
value of the LGL nodes at a particular element interface.
=#

# Point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# Reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for i in eachnode(semi.solver), j in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, j, element)
        # Changing the node positions passed to the initial condition by the minimum
        # amount possible with the current type of floating point numbers allows setting
        # discontinuous initial data in a simple way. In particular, a check like `if x < x_jump`
        # works if the jump location `x_jump` is at the position of an interface.
        if i == 1 && j == 1 # bottom left corner
            x_node = SVector(nextfloat(x_node[1]), nextfloat(x_node[2]))
        elseif i == 1 && j == nnodes(semi.solver) # top left corner
            x_node = SVector(nextfloat(x_node[1]), prevfloat(x_node[2]))
        elseif i == nnodes(semi.solver) && j == 1 # bottom right corner
            x_node = SVector(prevfloat(x_node[1]), nextfloat(x_node[2]))
        elseif i == nnodes(semi.solver) && j == nnodes(semi.solver) # top right corner
            x_node = SVector(prevfloat(x_node[1]), prevfloat(x_node[2]))
        elseif i == 1 # left boundary
            x_node = SVector(nextfloat(x_node[1]), x_node[2])
        elseif j == 1 # bottom boundary
            x_node = SVector(x_node[1], nextfloat(x_node[2]))
        elseif i == nnodes(semi.solver) # right boundary
            x_node = SVector(prevfloat(x_node[1]), x_node[2])
        elseif j == nnodes(semi.solver) # top boundary
            x_node = SVector(x_node[1], prevfloat(x_node[2]))
        end

        u_node = initial_condition_dam_break(x_node, first(tspan), equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
    end
end

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (energy_total,
                                                                 energy_kinetic,
                                                                 energy_internal))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 500,
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
summary_callback() # print the timer summary
