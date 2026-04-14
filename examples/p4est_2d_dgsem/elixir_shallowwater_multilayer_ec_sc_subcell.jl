using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations with a single layer and a bottom
# topography function for a blast wave test with discontinuous initial conditions to test entropy
# conservation with the subcell limiting volume integral on a curvilinear mesh.

equations = ShallowWaterMultiLayerEquations2D(gravity = 9.81, rhos = (1.0))

function initial_condition_weak_blast_wave(x, t,
                                           equations::ShallowWaterMultiLayerEquations2D)
    # Set up polar coordinates
    RealT = eltype(x)
    inicenter = SVector(convert(RealT, 0.0), convert(RealT, 0.0))
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    phi = atan(y_norm, x_norm)
    sin_phi, cos_phi = sincos(phi)

    # Calculate primitive variables
    H = r > 0.2f0 ? 3.8f0 : 4.0f0
    v1 = r > 0.2f0 ? zero(RealT) : convert(RealT, 0.1882) * cos_phi
    v2 = r > 0.2f0 ? zero(RealT) : convert(RealT, 0.1882) * sin_phi
    v1 = 0.0
    v2 = 0.0

    # Setup a continuous bottom topography
    r = 0.4
    b = (((x[1])^2 + (x[2])^2) < r^2 ?
         0.5 * (cos(1 / r * pi * sqrt((x[1])^2 + (x[2])^2)) + 1) : 0.0)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_weak_blast_wave

###############################################################################
# Get the DG approximation space

polydeg = 3
volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)

basis = LobattoLegendreBasis(polydeg)
# For reproducibility set the limiting coefficients to pure DG,
# see https://github.com/trixi-framework/Trixi.jl/pull/2007
limiter_idp = SubcellLimiterIDP(equations, basis;)
volume_integral = VolumeIntegralSubcellLimiting(limiter_idp;
                                                volume_flux_dg = volume_flux,
                                                volume_flux_fv = surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the P4estMesh and setup a periodic mesh on the domain [-1,1]^2 with an affine type mapping to
# obtain a warped curvilinear mesh.
function mapping_twist(xi, eta)
    y = eta + 0.1 * sin(pi * xi) * cos(0.5 * pi * eta)
    x = xi + 0.1 * sin(pi * eta) * cos(0.5 * pi * xi)
    return SVector(x, y)
end

# Create P4estMesh with 8 x 8 elements
trees_per_dimension = (2, 2)
mesh = P4estMesh(trees_per_dimension, polydeg = 2,
                 initial_refinement_level = 2,
                 periodicity = true,
                 mapping = mapping_twist)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver

tspan = (0.0, 0.2)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_polydeg = polydeg,
                                     extra_analysis_errors = (:conservation_error,))

stepsize_callback = StepsizeCallback(cfl = 0.5)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     extra_node_variables = (:limiting_coefficient,))

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################

# Setup the stage callbacks
stage_callbacks = (SubcellLimiterIDPCorrection(),)

# run the simulation
sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = stage_callbacks);
                  dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                  ode_default_options()...,
                  callback = callbacks, adaptive = false);
