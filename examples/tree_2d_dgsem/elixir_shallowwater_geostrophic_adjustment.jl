
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# semidiscretization of the shallow water equations with corriolis source term

equations = ShallowWaterEquationsWetDry2D(gravity_constant = 1.0, f0 = 2.0)

# Initial condition for a geostrophic adjustment test with elliptical mass imbalance. 
function initial_condition_geostrophic_adjustment(x, t, equations::ShallowWaterEquationsWetDry2D)
    # Parameters of the initial perturbation
    a_0 = 0.5    # amplitude
    r_e = 0.1   # edge width
    r_i = 1.0   # radius
    lambda = 2.5 # aspect ratio

    h = 1 + a_0 / 2 * (1 - tanh((sqrt((sqrt(lambda) * x[1])^2 + (x[2] / sqrt(lambda))^2) - r_i) / r_e))   
    h_v1 = 0.0
    h_v2 = 0.0
    b = 0.0

    return SVector(h, h_v1, h_v2, b)
end

# Outflow boundary condition for TreeMesh
@inline function boundary_condition_outflow(u_inner, orientation,
                                                    direction, x, t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterEquationsWetDry2D)

    h_inner = Trixi.waterheight(u_inner, equations)
    v1_inner, v2_inner = velocity(u_inner, equations)
    h_outer = 1.0

    surface_flux_function, _ = surface_flux_functions


    ## get the appropriate normal vector from the orientation
    if direction == 1
            v1_outer = v1_inner + 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
            u_boundary = SVector(h_outer, h_outer*v1_outer, u_inner[3], u_inner[4])
    elseif direction == 2
            v1_outer = v1_inner - 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
            u_boundary = SVector(h_outer, h_outer*v1_outer, u_inner[3], u_inner[4])
    elseif direction == 3
        v2_outer = v2_inner + 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
        u_boundary = SVector(h_outer, u_inner[2], h_outer*v2_outer, u_inner[4])
    elseif direction == 4
        v2_outer = v2_inner - 2 * sqrt(equations.gravity) * (sqrt(h_outer) - sqrt(h_inner))
        u_boundary = SVector(h_outer, u_inner[2], h_outer*v2_outer, u_inner[4])
    end

    walls = true

    if walls == true
        # Calculate boundary flux
        if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
            if (direction == 2 && sign(v1_inner)==1) || (direction == 4 && sign(v2_inner)==1)
                flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
            else
                flux, _ = boundary_condition_slip_wall(u_inner, orientation, direction, x, t, surface_flux_functions, equations)
            end
        else # u_boundary is "left" of boundary, u_inner is "right" of boundary
            if (direction == 1 && sign(v1_inner)==-1) || (direction == 3 && sign(v2_inner)==-1)
                flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
            else
            flux, _ = boundary_condition_slip_wall(u_inner, orientation, direction, x, t, surface_flux_functions, equations)
            end
        end
    else
        # Calculate boundary flux
        if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
                flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
        else # u_boundary is "left" of boundary, u_inner is "right" of boundary
                flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
        end
    end
    return (flux, zero(u_inner))
end

initial_condition = initial_condition_geostrophic_adjustment

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxPlusDissipation(flux_wintermeyer_etal, DissipationLocalLaxFriedrichs()), flux_nonconservative_wintermeyer_etal)

solver = DGSEM(polydeg = 3,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))


###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (-10.0, -10.0)
coordinates_max = (10.0, 10.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 40_000, periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, 
                                    source_terms = source_terms_corriolis, boundary_conditions=boundary_condition_outflow)

###############################################################################
# ODE solver

tspan = (0.0, 17.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_integrals = (lake_at_rest_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 1.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary
