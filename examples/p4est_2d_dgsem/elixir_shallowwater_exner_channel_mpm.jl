
using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Roots: Order2, solve, ZeroProblem

###############################################################################
# semidiscretization of the shallow water Exner equations with Meyer-Peter-Mueller
# sediment closure

equations = ShallowWaterExnerEquations2D(gravity = 9.81, H0 = 1.0,
                                         rho_f = 1000.0, rho_s = 2600.0, porosity = 0.4,
                                         friction = ManningFriction(n = 0.0196),
                                         sediment_model = MeyerPeterMueller(theta_c = 0.047,
                                                                            d_s = 1e-3))

# Sand dune in an L-shaped channel taken from Castro et al. https://doi.org/10.1016/j.cma.2009.03.001
@inline function initial_condition_L_shaped_channel(x, t,
                                                    equations::ShallowWaterExnerEquations2D)

    # Sediment bump
    h_b = 0.1 + 0.2 * exp(-((x[1] - 2)^2 + (x[2] - 2)^2))

    # Water height where H_0 = 1
    h = equations.H0 - h_b

    # Momenta
    q1 = zero(eltype(x))
    q2 = zero(eltype(x))

    return SVector(h, q1, q2, h_b)
end

initial_condition = initial_condition_L_shaped_channel

# From Castro et al. inflow needs v1=0, v2=-0.5, h_b=0.1 and outflow is "free".
# Here we set the inflow momentum to recover (approximately) these velocity values.
# At outflow, we set the water height. These boundary condition routines are modified
# versions of the standard shallow water case. Though not exactly characteristic BCs,
# they work because the sediment mode evolves very slowly and the shallow water part
# dominates.
function TrixiShallowWater.BoundaryConditionWaterHeight(h_boundary::Real,
                                                        equations::ShallowWaterExnerEquations2D{RealT}) where {RealT}
    # Convert function output to the correct type
    h_boundary = convert(RealT, h_boundary)
    return TrixiShallowWater.BoundaryConditionWaterHeight(t -> h_boundary)
end

function (boundary_condition::TrixiShallowWater.BoundaryConditionWaterHeight)(u_inner,
                                                                              normal_direction,
                                                                              x, t,
                                                                              surface_flux_functions,
                                                                              equations::ShallowWaterExnerEquations2D)
    surface_flux, nonconservative_flux_function = surface_flux_functions

    # Extract the gravitational acceleration
    g = equations.gravity

    # Normalized normal vector
    norm_ = Trixi.norm(normal_direction)
    normal = normal_direction / norm_

    # Apply the rotation that maps `normal` onto the x-axis to `u_inner`.
    u_rotated = Trixi.rotate_to_x(u_inner, normal, equations)

    # Get the water height and velocity from the inner state
    h_inner = waterheight(u_rotated, equations)
    v1, v2 = velocity(u_inner, equations)
    v_inner_normal, v_inner_tangential = velocity(u_rotated, equations)

    # Extract the external water height from the boundary condition and set sediment height
    h_boundary = boundary_condition.h_boundary(t)
    h_b = 0.1

    # In the case of inflow we fallback to setting a wall boundary condition.
    # TODO: Could try alternative of forcing outflow by taking the absolute value#
    # of the normal velocity and setting the external pressure like in FUN3D or FLEXI.
    # Big question is how to set the external pressure in this situation.
    if v_inner_normal < 0
        return Trixi.boundary_condition_slip_wall(u_inner, normal_direction, x, t,
                                                  surface_flux_functions, equations)
    else
        # Calculate the boundary state in the rotated coordinate system.
        # To extrapolate the external velocity assume that the Riemann invariant remains constant across
        # the incoming characteristic. In the case of inflow we assume that the tangential velocity at
        # the boundary is zero.
        v_boundary_normal = v_inner_normal - 2 * (sqrt(g * h_boundary) - sqrt(g * h_inner))
        v_boundary_normal < 0 ? hv_boundary_tangential = zero(u_rotated[3]) :
        hv_boundary_tangential = u_rotated[3]

        u_boundary = SVector(h_boundary, h_boundary * v_boundary_normal,
                             hv_boundary_tangential, h_b)

        # Compute the boundary flux in the rotated coordinate system using LLF
        local_lax_friedrichs_flux = FluxPlusDissipation(flux_ersing_etal,
                                                        DissipationLocalLaxFriedrichs())
        flux = local_lax_friedrichs_flux(u_rotated, u_boundary, 1, equations)
        noncons_flux = nonconservative_flux_function(u_rotated, u_boundary, 1, equations)

        # Apply the back-rotation that maps the x-axis onto `normal` to the boundary flux.
        flux = Trixi.rotate_from_x(flux, normal, equations) * norm_
        noncons_flux = Trixi.rotate_from_x(noncons_flux, normal, equations) * norm_

        # Return the conservative and nonconservative fluxes.
        return (flux, noncons_flux)
    end
end

function TrixiShallowWater.BoundaryConditionMomentum(hv1_boundary::Real, hv2_boundary::Real,
                                                     equations::ShallowWaterExnerEquations2D{RealT}) where {RealT}
    # Convert function output to the correct type
    hv1_boundary = convert(RealT, hv1_boundary)
    hv2_boundary = convert(RealT, hv2_boundary)
    return TrixiShallowWater.BoundaryConditionMomentum((t -> (hv1_boundary, hv2_boundary)))
end

function (boundary_condition::TrixiShallowWater.BoundaryConditionMomentum)(u_inner,
                                                                           normal_direction,
                                                                           x, t,
                                                                           surface_flux_functions,
                                                                           equations::ShallowWaterExnerEquations2D)
    _, nonconservative_flux_function = surface_flux_functions

    # Extract the gravitational acceleration
    g = equations.gravity

    # Normalized normal vector
    norm_ = Trixi.norm(normal_direction)
    normal = normal_direction / norm_

    # Apply the rotation that maps `normal` onto the x-axis to `u_inner` and `hv_boundary``.
    u_rotated = Trixi.rotate_to_x(u_inner, normal, equations)

    # Get the water height and velocity from the inner state
    h_inner = waterheight(u_rotated, equations)
    v_inner_normal, _ = velocity(u_rotated, equations)

    # Extract the external momentum from the boundary condition
    hv1_boundary, hv2_boundary = boundary_condition.hv_boundary(t)

    hv_boundary_normal = hv1_boundary * normal[1] + hv2_boundary * normal[2]
    hv_boundary_tangential = -hv1_boundary * normal[2] + hv2_boundary * normal[1]

    # Calculate the boundary state in the rotated coordinate system.
    # To extrapolate the external water height `h_boundary` assume that the Riemann invariant remains
    # constant across the incoming characteristic.
    # Requires one to solve for the roots of a nonlinear function, see Eq. (52) in the reference above.
    # For convenience we substitute x = h_boundary and solve for x.
    fx = ZeroProblem(x -> 2 * sqrt(g) * x^(3 / 2) -
                          (v_inner_normal + 2 * sqrt(g * h_inner)) * x +
                          hv_boundary_normal, h_inner)
    h_boundary = solve(fx, Order2())

    hv_boundary_normal < 0 ? nothing : hv_boundary_tangential = u_rotated[3]

    u_boundary = SVector(h_boundary, hv_boundary_normal, hv_boundary_tangential, u_inner[4])

    # Compute the boundary flux in the rotated coordinate system.
    flux = Trixi.flux(u_boundary, 1, equations)
    noncons_flux = nonconservative_flux_function(u_rotated, u_boundary, 1, equations)

    # Apply the back-rotation that maps the x-axis onto `normal` to the boundary flux.
    flux = Trixi.rotate_from_x(flux, normal, equations) * norm_
    noncons_flux = Trixi.rotate_from_x(noncons_flux, normal, equations) * norm_

    # Return the conservative and nonconservative fluxes.
    return (flux, noncons_flux)
end

boundary_condition_inflow = BoundaryConditionMomentum(0.0, -0.45, equations)
boundary_condition_outflow = BoundaryConditionWaterHeight(0.9, equations)

boundary_conditions = (; Channel = boundary_condition_slip_wall,
                       Inflow = boundary_condition_inflow,
                       Outflow = boundary_condition_outflow)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

surface_flux = (FluxPlusDissipation(flux_ersing_etal, dissipation_roe),
                flux_nonconservative_ersing_etal)

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = water_sediment_height)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the unstructured mesh from a file

default_meshfile = joinpath(@__DIR__, "mesh_L_channel.inp")

isfile(default_meshfile) ||
    Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/6026ce25bdf2a5498e43cbf3350c0397/raw/7bc7a519a728c07ba50ede815d7bbabc6313cbd0/mesh_L_channel.inp",
                   default_meshfile)

meshfile = default_meshfile

# Can run on a coarser mesh of 165 elements by setting `initial_refinement_level = 0`
mesh = P4estMesh{2}(meshfile, initial_refinement_level = 1)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_term_bottom_friction,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solver

tspan = (0.0, 3600.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 50.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks, maxiters = 1e9);
