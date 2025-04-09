
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break
# test over a dry domain with a discontinuous bottom topography function

equations = ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                              rhos = (0.8, 0.85, 0.9, 0.95, 1.0))

# Initial condition of a dam break with a discontinuous water heights and bottom topography.
# To test the wet/dry functionality this test case considers a dam break over a dry domain with
# multiple layers.
# Works as intended for TreeMesh1D with `initial_refinement_level=5`. If the mesh
# refinement level is changed the initial condition below may need changed as well to
# ensure that the discontinuities lie on an element interface.
function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations1D)
    # Set the discontinuity
    if x[1] <= 10.0
        H = [3.5, 3.0, 2.5, 2.0, 1.5]
        b = 0.0
    else
        # Right side of the domain is dry
        b = 0.5
        H = [b, b, b, b, b]
    end

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations1D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    for i in reverse(eachlayer(equations))
        if i == nlayers(equations)
            H[i] = max(H[i], b + equations.threshold_limiter)
        else
            H[i] = max(H[i], H[i + 1] + equations.threshold_limiter)
        end
    end

    # Initialize zero velocity
    v = zero(H)

    return prim2cons(SVector(H..., v..., b), equations)
end

initial_condition = initial_condition_dam_break

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(3)

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

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false)

stepsize_callback = StepsizeCallback(cfl = 0.2)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 500,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback, save_solution)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!);
            dt = 1.0, adaptive = false,
            ode_default_options()..., callback = callbacks);
