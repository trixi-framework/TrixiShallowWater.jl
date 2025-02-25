
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break test with a 
# discontinuous bottom topography function to test entropy conservation

equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 1.0,
                                              rhos = (0.9, 0.95, 1.0))

# This test case uses a special work around to setup a truly discontinuous bottom topography 
# function and initial condition for this academic testcase of entropy conservation. First, a 
# dummy initial_condition_dam_break is introduced to create the semidiscretization. Then the initial
# condition is reset with the true discontinuous values from initial_condition_discontinuous_dam_break.

function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 0.3 * exp(-0.5 * ((x[1] - sqrt(2) / 2)^2 + (x[2] - sqrt(2) / 2)^2))

    if x[1] < sqrt(2) / 2
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

boundary_condition_constant = BoundaryConditionDirichlet(initial_condition_dam_break)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
solver = DGSEM(polydeg = 6, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the unstructured quad mesh from a file (downloads the file if not available locally)
mesh_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/8f8cd23df27fcd494553f2a89f3c1ba4/raw/85e3c8d976bbe57ca3d559d653087b0889535295/mesh_alfven_wave_with_twist_and_flip.mesh",
                           joinpath(@__DIR__, "mesh_alfven_wave_with_twist_and_flip.mesh"))

mesh = UnstructuredMesh2D(mesh_file, periodicity = false)

# Boundary conditions
boundary_condition = Dict(:Top => boundary_condition_slip_wall,
                          :Left => boundary_condition_slip_wall,
                          :Right => boundary_condition_slip_wall,
                          :Bottom => boundary_condition_slip_wall)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
                                    solver, boundary_conditions = boundary_condition)

###############################################################################
# ODE solver

tspan = (0.0, 0.5)
ode = semidiscretize(semi, tspan)

###############################################################################
# Workaround to set a discontinuous bottom topography and initial condition for debugging and testing.

# alternative version of the initial conditinon used to setup a truly discontinuous
# test case and initial condition.
# In contrast to the usual signature of initial conditions, this one get passed the
# `element_id` explicitly. In particular, this initial conditions works as intended
# only for the specific mesh loaded above!
function initial_condition_discontinuous_dam_break(x, t, element_id,
                                                   equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 0.3 * exp(-0.5 * ((x[1] - sqrt(2) / 2)^2 + (x[2] - sqrt(2) / 2)^2))

    # Left side of discontinuity
    IDs = [1, 2, 5, 6, 9, 10, 13, 14]
    if element_id in IDs
        H = SVector(1.0, 0.8, 0.6)
    else # Right side of discontinuity
        H = SVector(0.9, 0.7, 0.5)
        b += 0.1
    end

    v1 = zero(H)
    v2 = zero(H)
    return prim2cons(SVector(H..., v1..., v2..., b),
                     equations)
end

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for j in eachnode(semi.solver), i in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, j, element)
        u_node = initial_condition_discontinuous_dam_break(x_node, first(tspan), element,
                                                           equations)
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
