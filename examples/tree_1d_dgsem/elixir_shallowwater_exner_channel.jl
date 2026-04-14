
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Roots

###############################################################################
# Semidiscretization of the shallow water exner equations for a channel flow problem
# with sediment transport

equations = ShallowWaterExnerEquations1D(gravity = 9.81, rho_f = 1.0, rho_s = 1.0,
                                         porosity = 0.4,
                                         sediment_model = GrassModel(A_g = 0.01))

# Initial condition for a channel flow problem over a sediment hump.
# An asymptotic solution based on the method of characteristics was derived under a rigid-lid
# approximation in chapter 3.5.1 of the thesis:
# - Justin Hudson (2001)
#   "Numerical Techniques for Morphodynamic Modelling"
#   [PhD Thesis, University of Reading](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=f78dbae9cfbb12ae823975c6ce9d2585b40417ba)
function initial_condition_channel(x, t,
                                   equations::ShallowWaterExnerEquations1D)
    (; sediment_model, porosity_inv) = equations

    H_ref = 10.0    # Reference water level
    hv = 10.0

    # Use the method of characteristics to find the asymptotic solution for the bed height, see
    # Eq. 3.16 in the reference.
    # First use an iterative method to predict x0
    f(x0) = x[1] - x0 -
            sediment_model.A_g * porosity_inv * 3 * hv^3 * t *
            (H_ref - sinpi((x0 - 300) / 200)^2)^(-4)
    fx = Roots.ZeroProblem(f, 400.0)
    x0 = Roots.solve(fx, atol = 0.0, rtol = 0.0)

    # If the result is outside 300 < x < 500 the result is invalid and instead compute x0 from
    if x0 > 500 || x0 < 300
        x0 = x[1] -
             sediment_model.A_g * porosity_inv * 3 * hv^3 * t * H_ref^(-4)
    end

    # Compute the sediment and water height
    300 < x0 < 500 ? h_b = 0.1 + sinpi((x0 - 300) / 200)^2 : h_b = 0.1
    h = H_ref - h_b

    return SVector(h, hv, h_b)
end

initial_condition = initial_condition_channel

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxPlusDissipation(flux_ersing_etal, dissipation_roe),
                flux_nonconservative_ersing_etal)

basis = LobattoLegendreBasis(6)

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
# Get the TreeMesh and setup a periodic mesh

coordinates_min = 0.0
coordinates_max = 1000.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver

tspan = (0.0, 30_000.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 0.8)

save_solution = SaveSolutionCallback(dt = 10_000.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback, save_solution)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
