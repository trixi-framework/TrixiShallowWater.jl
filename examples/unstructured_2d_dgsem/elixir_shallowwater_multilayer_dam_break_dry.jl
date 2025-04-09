
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break test over a dry domain
# with a discontinuous bottom topography function
equations = ShallowWaterMultiLayerEquations2D(gravity = 1.0,
                                              rhos = (0.9, 0.95, 1.0))

# This test case uses a special work around to setup a truly discontinuous bottom topography
# function and initial condition. First, a dummy initial_condition_dam_break is introduced to create
# the semidiscretization. Then the initial condition is reset with the true discontinuous values
# from initial_condition_discontinuous_dam_break.

function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 1.4 * exp(-10.0 * ((x[1] - sqrt(2) / 2)^2 + (x[2] - sqrt(2) / 2)^2))

    if x[1] < sqrt(2) / 2
        H = [1.0, 0.8, 0.6]
    else
        b += 0.1
        H = [b, b, b]
    end

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

    return prim2cons(SVector(H..., v1..., v2..., b),
                     equations)
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
basis = LobattoLegendreBasis(6)

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

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Workaround to set a discontinuous bottom topography and initial condition.

# alternative version of the initial conditinon used to setup a truly discontinuous
# test case and initial condition.
# In contrast to the usual signature of initial conditions, this one get passed the
# `element_id` explicitly. In particular, this initial conditions works as intended
# only for the specific mesh loaded above!
function initial_condition_discontinuous_dam_break(x, t, element_id,
                                                   equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 1.4 * exp(-10.0 * ((x[1] - sqrt(2) / 2)^2 + (x[2] - sqrt(2) / 2)^2))

    # Left side of discontinuity
    IDs = [1, 2, 5, 6, 9, 10, 13, 14]
    if element_id in IDs
        H = [1.0, 0.8, 0.6]
    else # Right side of discontinuity
        b += 0.1
        H = [b, b, b]
    end

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

analysis_interval = 10
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (energy_total,
                                                                 energy_kinetic,
                                                                 energy_internal))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);
