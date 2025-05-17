
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations with one layer
# to fallback to be the standard shallow water equations

equations = ShallowWaterMultiLayerEquations2D(gravity = 9.81, H0 = 0.0,
                                              threshold_limiter = 1e-12,
                                              rhos = (1.0))

"""
    initial_condition_three_mounds(x, t, equations::ShallowWaterMultiLayerEquations2D)

Initial condition simulating a dam break. The bottom topography is given by one large and two smaller
mounds. The mounds are flooded by the water for t > 0. To smooth the discontinuity, a logistic function
is applied.

The initial conditions is taken from Section 6.3 of the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and Timothy Warburton (2018)
  An entropy stable discontinuous Galerkin method for the shallow water equations on
  curvilinear meshes with wet/dry fronts accelerated by GPUs\n
  [DOI: 10.1016/j.jcp.2018.08.038](https://doi.org/10.1016/j.jcp.2018.08.038)
"""
function initial_condition_three_mounds(x, t, equations::ShallowWaterMultiLayerEquations2D)
    x1, x2 = x
    M_1 = 1 - 0.1 * sqrt((x1 - 30.0)^2 + (x2 - 22.5)^2)
    M_2 = 1 - 0.1 * sqrt((x1 - 30.0)^2 + (x2 - 7.5)^2)
    M_3 = 2.8 - 0.28 * sqrt((x1 - 47.5)^2 + (x2 - 15.0)^2)

    b = max(0.0, M_1, M_2, M_3)

    # use a logistic function to transfer water height value smoothly
    L = 1.875    # maximum of function
    x0 = 8  # center point of function
    k = -75.0 # sharpness of transfer

    H = [max(b, L / (1.0 + exp(-k * (x1 - x0))))]

    # Everything is initially not moving
    v1 = zero(H)
    v2 = zero(H)

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations2D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    for i in reverse(eachlayer(equations))
        if i == nlayers(equations)
            H[i] = max(H[i], b + equations.threshold_limiter)
        else
            H[i] = max(H[i], H[i + 1] + equations.threshold_limiter)
        end
    end

    # Return the conservative variables
    return prim2cons(SVector(H..., v1..., v2..., b), equations)
end

initial_condition = initial_condition_three_mounds

function boundary_condition_outflow(u_inner, normal_direction::AbstractVector, x, t,
                                    surface_flux_functions,
                                    equations::ShallowWaterMultiLayerEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions
    # Impulse and bottom from inside, height from external state
    u_outer = SVector(equations.threshold_limiter, u_inner[2], u_inner[3], u_inner[4])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

boundary_conditions = Dict(:Bottom => boundary_condition_slip_wall,
                           :Top => boundary_condition_slip_wall,
                           :Right => boundary_condition_outflow,
                           :Left => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(4)

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
# Get the unstructured quad mesh from a file (downloads the file if not available locally)

mesh_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/7a34a0116dc84b69449ea25054b7585a/raw/e098000588ac5d8b543fc807b7fcc49df63dc7a0/mesh_three_mound.inp",
                           joinpath(@__DIR__, "mesh_three_mound.inp"));

mesh = P4estMesh{2}(mesh_file)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solver

tspan = (0.0, 20.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,),
                                     save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.5,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.5)

# Cannot simply use `waterheight` here for multilayer equations.
# Need a helper function to extract the relevant variable.
# TODO: This AMR indicator fires too often in the dry regions
@inline function main_waterheight(u, equations)
    return waterheight(u, equations)[1]
end

amr_indicator = IndicatorLöhner(semi, variable = main_waterheight)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 0,
                                      med_level = 1, med_threshold = 0.1,
                                      max_level = 3, max_threshold = 0.25)

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 5,
                           adapt_initial_condition = false,
                           adapt_initial_condition_only_refine = false)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        amr_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()...,
            callback = callbacks,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            adaptive = false);
