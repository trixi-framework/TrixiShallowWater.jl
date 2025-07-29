
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Roots

###############################################################################
# semidiscretization of the shallow water equations for a transonic moving water 
# steady-state with a standing shock.

equations = ShallowWaterEquations1D(gravity = 9.812, H0 = 3.25)

"""
    inverse_transform(E, hv, sigma, b)

Inverse transformation from equilibrium variables (E, hv) to conservative variables (h, hv). Besides the
equilibrium variables, which are the total energy `E` and momentum `hv`, the function also depends 
on the bottom topography `b` and the flow regime `sigma` (supersonic = 1 , sonic = 0 or subsonic = -1).

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

    # Check if the state is admissible and compute the water height
    # as in equation (2.19) of the reference in the docstring.
    if sigma == 0 && E_hat ≈ 1 # sonic state
        h_hat = 1
    elseif abs(sigma) == 1 && E_hat > 1 # supersonic / subsonic state
        if sigma == 1 # supersonic
            h_hat_init = 0.5 # needs to be < 1
        else # subsonic
            h_hat_init = 2.0 # needs to be > 1
        end

        # Setup the root finding problem.
        f(h_hat) = E_hat - 2 / 3 * ((1 / (2 * h_hat^2)) + h_hat)
        D(f) = h_hat -> Trixi.ForwardDiff.derivative(f, float(h_hat))

        # Solve the root finding problem using Newton's method
        h_hat = Roots.newton((f, D(f)), h_hat_init)
    else
        throw(error("The given state is not admissible: E_hat = $E_hat, sigma = $sigma"))
    end

    # Compute and return the water height `h` from the normalized water height `h_hat = h / h_0`.
    return h_hat * h_0
end

"""
    initial_condition_moving_water_transonic(x, t, equations::ShallowWaterEquations1D)

Set the initial condition for a transonic moving water steady-state and a quadratic bottom 
topography, to test the well-balancedness of the scheme.

The test parameters are taken from Section 4.1 of the paper:
- Sebastian Noelle, Yulong Xing and Chi-Wang Shu (2007)
  High Order Well-balanced Finite Volume WENO Schemes for Shallow Water Equation with Moving Water
  [DOI: 10.1016/j.jcp.2007.03.031](https://doi.org/10.1016/j.jcp.2007.03.031).
"""
function initial_condition_moving_water_shock(x, t,
                                              equations::ShallowWaterEquations1D)
    # Extract the gravitational acceleration
    g = equations.gravity

    hv = 0.18 # momentum

    # Set the total energy before and after the shock
    if x[1] < 11.665504281554291
        E = 3 / 2 * (g * 0.18)^(2 / 3) + g * 0.2
    else
        E = 0.18^2 / (2 * 0.33^2) + g * 0.33
    end

    # # Set the quadratic bottom topography function
    if 8 <= x[1] <= 12
        b = 0.2 - 0.05 * (x[1] - 10.0)^2
    else
        b = 0.0
    end

    # Set the sign function to label the flow regime (subsonic = -1, sonic = 0, supersonic = 1).
    # A small tolerance is required to avoid numerical issues in the inverse_transform function
    # close to the sonic point at x = 10.
    tol = 1e-12
    if x[1] <= 10.0 - tol || x[1] >= 11.665504281554291
        sigma = -1 # subsonic
    elseif 10 - tol < x[1] < 10.0 + tol
        sigma = 0 # sonic
    elseif 10 + tol <= x[1] < 11.665504281554291
        sigma = 1 # supersonic
    end

    h = inverse_transform(E, hv, sigma, b)

    return SVector(h, hv, b)
end

initial_condition = initial_condition_moving_water_shock

boundary_condition_inflow = BoundaryConditionMomentum(0.18, equations)
boundary_condition_outflow = BoundaryConditionDirichlet(initial_condition_moving_water_shock)

boundary_conditions = (x_neg = boundary_condition_inflow,
                       x_pos = boundary_condition_outflow)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
# Up to Trixi.jl version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `DissipationLocalLaxFriedrichs(), i.e., `DissipationLocalLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 of Trixi.jl the default wave speed used in the LLF flux and dissipation operator to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `DissipationLocalLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.
surface_flux = (FluxPlusDissipation(flux_wintermeyer_etal,
                                    DissipationLocalLaxFriedrichs(max_abs_speed_naive)),
                flux_nonconservative_wintermeyer_etal)

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the TreeMesh and setup a non-periodic mesh

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

stepsize_callback = StepsizeCallback(cfl = 0.7)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary
