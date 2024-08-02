using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the SWE-Exner equations with source terms for convergence testing

# Equations with MeyerPeterMueller model
equations_mpm = ShallowWaterExnerEquations1D(gravity_constant = 10.0, rho_f = 0.5,
                                             rho_s = 1.0, porosity = 0.5,
                                             friction = ManningFriction(n = 0.01),
                                             sediment_model = MeyerPeterMueller(theta_c = 0.0,
                                                                                d_s = 1e-3))

# Equations with Grass model
equations_grass = ShallowWaterExnerEquations1D(gravity_constant = 10.0, rho_f = 0.5,
                                               rho_s = 1.0, porosity = 0.5,
                                               friction = ManningFriction(n = 0.0),
                                               sediment_model = GrassModel(A_g = 0.01))

equations = equations_grass

# Smooth initial condition to test convergence
@inline function Trixi.initial_condition_convergence_test(x, t,
                                                          equations::ShallowWaterExnerEquations1D)
    ω = sqrt(2) * pi

    h = 2.0 + cos(ω * x[1]) * cos(ω * t)
    v = 0.5
    h_b = 2.0 + sin(ω * x[1]) * cos(ω * t)

    return SVector(h, h * v, h_b)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterExnerEquations1D{T, S, ShieldsStressModel{T}}) where {T, S})

Source terms used for convergence tests in combination with [`Trixi.initial_condition_convergence_test`](@extref) 
when using the [`MeyerPeterMueller`](@ref) model.
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterExnerEquations1D{T,
                                                                                             S,
                                                                                             ShieldsStressModel{T}}) where {
                                                                                                                            T,
                                                                                                                            S
                                                                                                                            }
    ω = sqrt(2.0) * pi
    (; gravity, porosity_inv, rho_f, rho_s, r) = equations

    n = equations.friction.n

    # Constant expression from the MPM model
    c = sqrt(gravity * (1 / r - 1)) * 8.0 * porosity_inv *
        (rho_f / (rho_s - rho_f))^(3 / 2) * n^3

    h = -cos(x[1] * ω) * sin(t * ω) * ω - 0.5 * sin(x[1] * ω) * cos(t * ω) * ω

    hv = ((5.0 * c *
           (cos(x[1] * ω) * cos(t * ω) * ω - 0.5 * sin(x[1] * ω) * cos(t * ω) * ω)) /
          ((2.0 + cos(x[1] * ω) * cos(t * ω))^0.5) - 0.5 * cos(x[1] * ω) * sin(t * ω) * ω -
          0.25 * sin(x[1] * ω) * cos(t * ω) * ω +
          10.0 * (2.0 + cos(x[1] * ω) * cos(t * ω)) *
          (cos(x[1] * ω) * cos(t * ω) * ω - sin(x[1] * ω) * cos(t * ω) * ω))

    h_b = ((0.5 * ((0.125 * c) / (2.0 + cos(x[1] * ω) * cos(t * ω))) * sin(x[1] * ω) *
            cos(t * ω) * ω) / ((2.0 + cos(x[1] * ω) * cos(t * ω))^0.5) -
           sin(x[1] * ω) * sin(t * ω) * ω)

    return SVector(h, hv, h_b)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterExnerEquations1D{T, S, GrassModel{T}}) where {T, S})

Source terms used for convergence tests in combination with [`initial_condition_convergence_test`](@extref) 
when using the the [`GrassModel`](@ref) model.
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterExnerEquations1D{T,
                                                                                             S,
                                                                                             GrassModel{T}}) where {
                                                                                                                    T,
                                                                                                                    S
                                                                                                                    }
    ω = sqrt(2.0) * pi
    A_g = equations.sediment_model.A_g

    h = -cos(x[1] * ω) * sin(t * ω) * ω - 0.5 * sin(x[1] * ω) * cos(t * ω) * ω
    hv = -0.5 * cos(x[1] * ω) * sin(t * ω) * ω - 0.25 * sin(x[1] * ω) * cos(t * ω) * ω +
         10.0 * A_g *
         (cos(x[1] * ω) * cos(t * ω) * ω - 0.5 * sin(x[1] * ω) * cos(t * ω) * ω) +
         10.0 * (2.0 + cos(x[1] * ω) * cos(t * ω)) *
         (cos(x[1] * ω) * cos(t * ω) * ω - sin(x[1] * ω) * cos(t * ω) * ω)
    h_b = -sin(x[1] * ω) * sin(t * ω) * ω
    return SVector(h, hv, h_b)
end

initial_condition = initial_condition_convergence_test

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

solver = DGSEM(polydeg = 4,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = 0.0
coordinates_max = sqrt(2.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 2,
                n_cells_max = 10_000,
                periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_test)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 200
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

save_solution = SaveSolutionCallback(interval = 200,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback, save_solution)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);

summary_callback() # print the timer summary
