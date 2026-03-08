using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Symbolics

###############################################################################
# Semidiscretization of the multilayer shallow water equations with a single layer to test
# convergence. The initial condition and source term are created using the 
# method of manufactured solutions (MMS) with the help of Symbolics.jl.

equations = ShallowWaterMultiLayerEquations2D(gravity = 1.1, rhos = (1.0))

### Create manufactured solution for method of manufactured solutions (MMS)

# Symbolic Variables
@variables x_sym[1:2], t_sym, g

# Define Differentials
Dt, Dx, Dy = Differential(t_sym), Differential(x_sym[1]), Differential(x_sym[2])

## Initial condition
###############################################################################
# Primitive Variables
ω = pi * 1.0
H = 4.0 + 0.2 * cos(ω * x_sym[1] + t_sym) + 0.2 * cos(ω * x_sym[2] + t_sym)
v1 = 0.5
v2 = 0.5
b = 1.0 + 0.2 * cos(ω * x_sym[1]) + 0.2 * cos(ω * x_sym[2])
h = H - b

init = [H, v1, v2, b]

## PDE
###############################################################################
eqs = [
    Dt(h) + Dx(h * v1) + Dy(h * v2),
    Dt(h * v1) + Dx(h * v1^2 + 0.5 * g * h^2) + Dy(h * v1 * v2) + g * h * Dx(b),
    Dt(h * v2) + Dx(h * v1 * v2) + Dy(h * v2^2 + 0.5 * g * h^2) + g * h * Dy(b),
    0
]

## Create the functions for the manufactured solution
########################################################################################
# Expand derivatives
du_exprs = expand_derivatives.(eqs)

# Build functions
du_funcs = build_function.(du_exprs, Ref(x_sym), t_sym, g, expression = Val(false))

init_funcs = build_function.(init, Ref(x_sym), t_sym, expression = Val(false))

function initial_condition_convergence_mms(x,
                                           t,
                                           equations::ShallowWaterMultiLayerEquations2D)
    prim = SVector{4, Float64}([f(x, t) for f in init_funcs]...)
    return prim2cons(prim, equations)
end

function source_terms_convergence_mms(u,
                                      x,
                                      t,
                                      equations::ShallowWaterMultiLayerEquations2D)
    g = equations.gravity
    return SVector{4, Float64}([f(x, t, g) for f in du_funcs]...)
end

initial_condition = initial_condition_convergence_mms

###############################################################################
# Get the DG approximation space

polydeg = 3
volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal_local_jump)

basis = LobattoLegendreBasis(polydeg)
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

# Create warped P4estMesh with 8 x 8 elements
trees_per_dimension = (2, 2)
mesh = P4estMesh(trees_per_dimension, polydeg = 3,
                 initial_refinement_level = 2,
                 periodicity = true,
                 mapping = mapping_twist)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence_mms,
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.1)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 500,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.7)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

stage_callbacks = (SubcellLimiterIDPCorrection(),)

sol = Trixi.solve(ode, Trixi.SimpleSSPRK33(stage_callbacks = stage_callbacks);
                  dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                  ode_default_options()...,
                  callback = callbacks, adaptive = false);
