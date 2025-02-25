
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break
# test with a discontinuous bottom topography function to test entropy conservation

equations = ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                              rhos = (0.85, 0.9, 1.0))

# Initial condition of a dam break with a discontinuous water heights and bottom topography.
# Works as intended for TreeMesh1D with `initial_refinement_level=5`. If the mesh
# refinement level is changed the initial condition below may need changed as well to
# ensure that the discontinuities lie on an element interface.
function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations1D)
    v = [0.0, 0.0, 0.0]

    # Set the discontinuity
    if x[1] <= 10.0
        H = [6.0, 4.0, 2.0]
        b = 0.0
    else
        H = [5.5, 3.5, 1.5]
        b = 0.5
    end

    return prim2cons(SVector(H..., v..., b), equations)
end

initial_condition = initial_condition_dam_break

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
solver = DGSEM(polydeg = 3,
               surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a non-periodic mesh

coordinates_min = 0.0
coordinates_max = 20.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 10000,
                periodicity = false)

boundary_condition = boundary_condition_slip_wall

# create the semidiscretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers

tspan = (0.0, 0.4)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (energy_total,
                                                                 energy_kinetic,
                                                                 energy_internal))

stepsize_callback = StepsizeCallback(cfl = 1.0)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 500,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution)

###############################################################################
# run the simulation

# use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL49(); abstol = 1.0e-8, reltol = 1.0e-8,
            ode_default_options()..., callback = callbacks);
summary_callback() # print the timer summary
