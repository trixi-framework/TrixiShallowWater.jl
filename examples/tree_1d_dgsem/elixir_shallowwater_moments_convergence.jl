using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Symbolics

###############################################################################
# Semidiscretization of the shallow water moment equations to test convergence against a 
# manufactured solution.

equations = ShallowWaterMomentEquations1D(gravity = 9.81, n_moments = 2)

### Create manufactured solution for method of manufactured solutions (MMS)
n = equations.n_moments

# Symbolic Variables
@variables x_sym, t_sym, g

# Define Differentials
Dt, Dx = Differential(t_sym), Differential(x_sym)

##  Initial condition
##################################################################################################
H = 7 + cos(sqrt(2) * 2 * pi * x_sym) * cos(2 * pi * t_sym)
v = 0.5
b = 2 + 0.5 * sinpi(sqrt(2) * x_sym)
h = H - b
a = [0.5 for i in 1:n]

init = [H, v, a..., b]

##  PDE 
###################################################################################################
# precompute the sum term
sum_moments = sum(h * a[j]^2 / (2j + 1) for j in 1:n)
sum_A = [sum(h * equations.A[i, j, k] * a[j] * a[k] for j in 1:n, k in 1:n) for i in 1:n]
sum_B = [sum(equations.B[i, j, k] * a[k] * Dx(h * a[j]) for j in 1:n, k in 1:n)
         for i in 1:n]

# additional moment equations 
mom_eqs = [Dt(h * a[i]) + Dx(2 * h * v * a[i] + sum_A[i]) - v * Dx(h * a[i]) + sum_B[i]
           for i in 1:n]

# PDE Source Terms
eqs = [
    Dt(h) + Dx(h * v),
    Dt(h * v) + Dx(h * v^2 + sum_moments) + g * h * Dx(h + b),
    mom_eqs...,
    0
]

## Create the functions for the manufactured solution
###################################################################################################
# Expand derivatives
du_exprs = expand_derivatives.(eqs)

# Build functions
du_funcs = build_function.(du_exprs, Ref(x_sym), Ref(t_sym), g, expression = Val(false))

init_funcs = build_function.(init, Ref(x_sym), t_sym, expression = Val(false))

# Trixi functions
function initial_condition_convergence(x, t, equations::ShallowWaterMomentEquations1D)
    prim = SVector{3 + n, Float64}([f(x[1], t) for f in init_funcs]...)
    return prim2cons(prim, equations)
end

function source_terms_convergence(u, x, t, equations::ShallowWaterMomentEquations1D)
    g = equations.gravity
    return SVector{3 + n, Float64}([f(x[1], t, g) for f in du_funcs]...)
end

initial_condition = initial_condition_convergence

###############################################################################
# Get the DG approximation space
volume_flux = (flux_careaga_etal, flux_nonconservative_careaga_etal)
surface_flux = (FluxPlusDissipation(flux_careaga_etal,
                                    DissipationLaxFriedrichsEntropyVariables(max_abs_speed)),
                flux_nonconservative_careaga_etal)

solver = DGSEM(polydeg = 3,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = 0.0
coordinates_max = sqrt(2.0)
mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000,
                periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh,
                                    equations,
                                    initial_condition,
                                    solver,
                                    source_terms = source_terms_convergence;
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.05)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 200,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution)

###############################################################################
# run the simulation

# use a Runge-Kutta method fixed step size
sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false);
            dt = 1e-4,
            ode_default_options()...,
            callback = callbacks,);
