
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break test over a dry domain
# with a discontinuous bottom topography function.
equations = ShallowWaterMultiLayerEquations2D(gravity = 1.0,
                                              rhos = (0.9, 0.95, 1.0))

function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 1.4 * exp(-10.0 * (x[1]^2 + x[2]^2))
    if x[1] > 0.0
        b += 0.1
    end

    if x[1] < -0.5
        H = [1.0, 0.8, 0.6]
    else
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

boundary_conditions = boundary_condition_slip_wall

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
# Get the TreeMesh and setup a non-periodic mesh with wall boundary conditions

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
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

analysis_interval = 50
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (energy_total,
                                                                 energy_kinetic,
                                                                 energy_internal))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation

# use a Runge-Kutta method with CFL-based time step
sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);
