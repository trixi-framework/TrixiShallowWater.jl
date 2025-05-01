using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Roots

###############################################################################
# semidiscretization of the shallow water equations for a subsonic moving water steady-state

equations = ShallowWaterEquationsWetDry1D(gravity = 9.812, H0 = 3.25)

"""
    inverse_transform(E, hv, sigma, b)

Inverse transformation from equilibrium variables (E,hv) to conservative variables (h,hv). Besides the
equilibrium variables, which are the total energy `E` and momentum `hv`, the function also depends 
on the bottom topography `b` and the flow regime `sigma`(supersonic = 1 , sonic = 0 or subsonic = -1).

The implementation follows the procedure described in Section 2.1 of the paper:
- Sebastian Noelle, Yulong Xing and Chi-Wang Shu (2007)
    High Order Well-balanced Finite Volume WENO Schemes for Shallow Water Equation with Moving Water
    [DOI: 10.1016/j.jcp.2007.03.031](https://doi.org/10.1016/j.jcp.2007.03.031).
"""
function inverse_transform(E, hv, sigma, b)
    # Extract the gravitational acceleration
    g = equations.gravity

    # Compute water height and specific energy at the sonic point
    h_0 = 1 / g * (g * abs(hv))^(2 / 3)
    phi_0 = 3 / 2 * (g * abs(hv))^(2 / 3)

    # normalized total energy
    E_hat = (E - g * b) / phi_0

    # Check if the state is admissible and compute the water height from equation (2.9).
    if (sigma == 0) && (E_hat ≈ 1) # sonic state
        h_hat = 1
    elseif abs(sigma) == 1 && E_hat > 1 # supersonic / subsonic state
        # Pick an initial guess for the root finding problem based on the flow regime
        if sigma == 1   # supersonic
            h_hat_init = 0.5 # needs to be < 1
        else    # subsonic
            h_hat_init = 2 # needs to be > 1
        end

        # Setup the root finding problem.
        f(h_hat) = E_hat - 2 / 3 * ((1 / (2 * h_hat^2)) + h_hat)
        D(f) = h_hat -> Trixi.ForwardDiff.derivative(f, float(h_hat))

        # Solve the root finding problem using Newton's method
        h_hat = Roots.newton((f, D(f)), h_hat_init)
    else
        throw(error("The given state is not admissible: E_hat = $E_hat, sigma = $sigma"))
    end

    # Return the water height `h = h_hat * h_0`.
    return h_hat * h_0
end

"""
    initial_condition_moving_water_subsonic(x, t, equations::ShallowWaterEquations1D)

Set the intitial condition for a subsonic moving water steady-state and quadratic bottom topography,
to test the well-balancedness of the scheme.

The test setup is taken from Section 4.1 of the paper:
- Sebastian Noelle, Yulong Xing and Chi-Wang Shu (2007)
    High Order Well-balanced Finite Volume WENO Schemes for Shallow Water Equation with Moving Water
    [DOI: 10.1016/j.jcp.2007.03.031](https://doi.org/10.1016/j.jcp.2007.03.031).
"""
function initial_condition_moving_water_subsonic(x, t,
                                                 equations::ShallowWaterEquationsWetDry1D)
    # Set initial conditions
    hv = 4.42 # momentum
    E = 22.06605 # total energy

    # Set the quadratic bottom topography function
    if 8 <= x[1] <= 12
        b = 0.2 - 0.05(x[1] - 10.0)^2
    else
        b = 0.0
    end

    sigma = -1 # sign function to label the flow regime (subsonic = -1, sonic = 0, supersonic = 1)

    # Compute the water height using the inverse transformation
    h = inverse_transform(E, hv, sigma, b)

    return SVector(h, hv, b)
end

initial_condition = initial_condition_moving_water_subsonic

boundary_condition_inflow = BoundaryConditionMomentum(4.42, equations)
boundary_condition_outflow = BoundaryConditionWaterHeight(2.0, equations)

boundary_conditions = (x_neg = boundary_condition_inflow,
                       x_pos = boundary_condition_outflow)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)

solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = 0.0
coordinates_max = 32.0  # This needs to be a multiple of 2 to match the corners of the bottom topography
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 10_000,
                periodicity = false)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solver

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = true)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary
