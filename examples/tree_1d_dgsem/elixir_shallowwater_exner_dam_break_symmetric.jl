
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the shallow water exner equations for a subcritical symmetric dam break problem

equations = ShallowWaterExnerEquations1D(gravity = 9.81, rho_f = 0.3, rho_s = 1.0,
                                         porosity = 0.4,
                                         sediment_model = GrassModel(A_g = 0.01))

# Initial condition for a synthetic test case of a symmetric dam break problem.
# The setup is taken from the paper:
# - S. Martínez-Aranda, J. Murillo, P. García-Navarro (2021)
#   "Comparison of new efficient 2D models for the simulation of bedload transport using the augmented
#    Roe approach"
#   [DOI: 10.1016/j.advwatres.2021.103931](https://doi.org/10.1016/j.advwatres.2021.103931)
function initial_condition_dam_break_symmetric(x, t,
                                               equations::ShallowWaterExnerEquations1D)
    # Setup initial perturbation of the water height
    if -0.5 <= x[1] <= 0.5
        h = 1.0
    else
        h = 0.2
    end

    # Set constant values for the sediment height and zero momentum
    hv = 0.0
    h_b = 1.0

    return SVector(h, hv, h_b)
end

initial_condition = initial_condition_dam_break_symmetric

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
# Get the TreeMesh and setup a periodic mesh

coordinates_min = -4.0
coordinates_max = 4.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solver

tspan = (0.0, 0.6)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback, save_solution)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
