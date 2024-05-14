
using Downloads: download
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations with a continuous
# bottom topography function

equations = ShallowWaterEquationsWetDry2D(gravity_constant = 9.812, H0 = 3.0)

function initial_condition_wb_testing(x, t, equations::ShallowWaterEquationsWetDry2D)
    # Set up polar coordinates
    inicenter = SVector(0.15, 0.15)
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)

    # Calculate primitive variables
    v1 = 0.0
    v2 = 0.0
    # v1 = r < 0.6 ? 1.75 : 0.0
    # v2 = r < 0.6 ? -1.75 : 0.0
    # bottom topography taken from Pond.control in [HOHQMesh](https://github.com/trixi-framework/HOHQMesh)
    x1, x2 = x
    b = (1.5 / exp(0.5 * ((x1 - 1.0)^2 + (x2 - 1.0)^2))
         +
         0.75 / exp(0.5 * ((x1 + 1.0)^2 + (x2 + 1.0)^2)))
    # H = equations.H0
    H = max(equations.H0, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_wb_testing

boundary_condition = Dict(:OuterCircle => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

# # ersing flux in both
# volume_flux = (flux_wintermeyer_etal, flux_nonconservative_ersing_etal)
# surface_fluxes = (flux_wintermeyer_etal, flux_nonconservative_ersing_etal)

# Wintermeyer flux in the volume
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

# # Wintermeyer flux in the surface for testing
# surface_fluxes = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

surface_fluxes = (FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
                  flux_nonconservative_chen_noelle)


# # Fjordholms flux for testing
# surface_fluxes = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal)

# Audusse HR with Rusanov flux
# surface_fluxes = (FluxHydrostaticReconstruction(flux_lax_friedrichs,
#                                                 hydrostatic_reconstruction_audusse_etal),
#                   flux_nonconservative_audusse_etal)

# Audusse HR with HLL flux
# surface_fluxes = (FluxHydrostaticReconstruction(flux_hll, hydrostatic_reconstruction_audusse_etal),
#                   flux_nonconservative_audusse_etal)

# Create the solver
solver = DGSEM(polydeg = 4, surface_flux = surface_fluxes,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# This setup is for the curved, split form well-balancedness testing

# Get the unstructured quad mesh from a file (downloads the file if not available locally)

default_mesh_file = joinpath(@__DIR__, "abaqus_outer_circle.inp")
isfile(default_mesh_file) ||
    download("https://gist.githubusercontent.com/andrewwinters5000/df92dd4986909927e96af23c37f6db5f/raw/8620823342f98c505a36351b210aab7f3b368041/abaqus_outer_circle.inp",
             default_mesh_file)
mesh_file = default_mesh_file

mesh = P4estMesh{2}(mesh_file)

# Refine bottom left quadrant of each tree to level 2
function refine_fn(p4est, which_tree, quadrant)
    quadrant_obj = unsafe_load(quadrant)
    if quadrant_obj.x == 0 && quadrant_obj.y == 0 && quadrant_obj.level < 2
          # return true (refine)
        return Cint(1)
    else
        # return false (don't refine)
        return Cint(0)
    end
end

# Refine recursively until each bottom left quadrant of a tree has level 2
# The mesh will be rebalanced before the simulation starts
refine_fn_c = @cfunction(refine_fn, Cint,
                         (Ptr{Trixi.p4est_t}, Ptr{Trixi.p4est_topidx_t},
                          Ptr{Trixi.p4est_quadrant_t}))
Trixi.refine_p4est!(mesh.p4est, true, refine_fn_c, C_NULL)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks, etc.

tspan = (0.0, 3.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# # Workaround to set a discontinuous bottom topography for debugging and testing.

# # alternative version of the initial conditinon used to setup a truly discontinuous
# # bottom topography function for this academic testcase.
# # The errors from the analysis callback are not important but the error for this lake at rest test case
# # `∑|H0-(h+b)|` should be around machine roundoff
# # In contrast to the usual signature of initial conditions, this one get passed the
# # `element_id` explicitly. In particular, this initial conditions works as intended
# # only for the specific mesh loaded above!
# function initial_condition_discontinuous_well_balancedness(x, t, element_id, equations::ShallowWaterEquations2D)
#   # Set the background values
#   H = equations.H0
#   v1 = 0.0
#   v2 = 0.0

#   x1, x2 = x
#   b = (  1.5 / exp( 0.5 * ((x1 - 1.0)^2 + (x2 - 1.0)^2) )
#        + 0.75 / exp( 0.5 * ((x1 + 1.0)^2 + (x2 + 1.0)^2) ) )

#   # Setup a discontinuous bottom topography using the element id number
# #  if element_id > 50 && element_id < 80 # for the original mesh file
# #  if element_id > 200 && element_id < 300 # for no hanging node grid with < 1 in function above
#   if element_id > 300 && element_id < 400 # for the forced hanging node grid with < 2 in function above
# #  if element_id > 400 && element_id < 650 # for the forced hanging node grid with < 3 in function above
#     b = 0.0
#   end

#   return prim2cons(SVector(H, v1, v2, b), equations)
# end

# # point to the data we want to augment
# u = Trixi.wrap_array(ode.u0, semi)
# # reset the initial condition
# for element in eachelement(semi.solver, semi.cache)
#   for j in eachnode(semi.solver), i in eachnode(semi.solver)
#     x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations, semi.solver, i, j, element)
#     u_node = initial_condition_discontinuous_well_balancedness(x_node, first(tspan), element, equations)
#     Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
#   end
# end
###############################################################################

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 1.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary
