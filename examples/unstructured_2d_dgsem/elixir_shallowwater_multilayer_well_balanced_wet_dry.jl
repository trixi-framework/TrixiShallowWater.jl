
using Downloads: download
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations with a discontinuous
# bottom topography function (set in the initial conditions) and dry domains

equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81, H0 = 1.5,
                                              rhos = (0.9, 0.95, 1.0))

# An initial condition with constant total water height and zero velocities to test well-balancedness.
# Note, this routine is used to compute errors in the analysis callback but the initialization is
# overwritten by `initial_condition_discontinuous_well_balancedness` below.
function initial_condition_well_balancedness(x, t,
                                             equations::ShallowWaterMultiLayerEquations2D)
    # Set the background values
    H = [1.5, 1.0, 0.5]
    v1 = zero(H)
    v2 = zero(H)

    # bottom topography taken from Pond.control in [HOHQMesh](https://github.com/trixi-framework/HOHQMesh)
    x1, x2 = x
    b = (1.5 / exp(0.5 * ((x1 - 1.0)^2 + (x2 - 1.0)^2)) +
         0.75 / exp(0.5 * ((x1 + 1.0)^2 + (x2 + 1.0)^2)))

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

    return prim2cons(SVector(H..., v1..., v2..., b), equations)
end

initial_condition = initial_condition_well_balancedness

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
# This setup is for the curved, split form well-balancedness testing

# Get the unstructured quad mesh from a file (downloads the file if not available locally)
default_mesh_file = joinpath(@__DIR__, "mesh_alfven_wave_with_twist_and_flip.mesh")
isfile(default_mesh_file) ||
    download("https://gist.githubusercontent.com/andrewwinters5000/8f8cd23df27fcd494553f2a89f3c1ba4/raw/85e3c8d976bbe57ca3d559d653087b0889535295/mesh_alfven_wave_with_twist_and_flip.mesh",
             default_mesh_file)
mesh_file = default_mesh_file

mesh = UnstructuredMesh2D(mesh_file, periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solver

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Workaround to set a discontinuous bottom topography for debugging and testing.

# alternative version of the initial conditinon used to setup a truly discontinuous
# bottom topography function for this academic testcase.
# The errors from the analysis callback are not important but the error for this lake at rest test case
# `∑|H0-(h+b)|` should be around machine roundoff
# In contrast to the usual signature of initial conditions, this one get passed the
# `element_id` explicitly. In particular, this initial conditions works as intended
# only for the specific mesh loaded above!
function initial_condition_discontinuous_well_balancedness(x, t, element_id,
                                                           equations::ShallowWaterMultiLayerEquations2D)
    # Set the background values
    H = [1.5, 1.0, 0.5]
    v1 = zero(H)
    v2 = zero(H)

    b = 0.2 + 0.1 * sin(2.0 * pi * x[1]) + 0.1 * cos(2.0 * pi * x[2])

    # Setup a discontinuous bottom topography to create wet/dry transitions using the `element_id` number
    if element_id == 7
        b += 0.1
    elseif element_id == 6
        b += 0.5
    elseif element_id == 10
        b += 1.0
    elseif element_id == 11
        b += 1.5
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

    return prim2cons(SVector(H..., v1..., v2..., b), equations)
end

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for j in eachnode(semi.solver), i in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, j, element)
        u_node = initial_condition_discontinuous_well_balancedness(x_node, first(tspan),
                                                                   element, equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
    end
end

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);

summary_callback() # print the timer summary
