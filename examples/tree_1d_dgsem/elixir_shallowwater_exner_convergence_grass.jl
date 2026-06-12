using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Symbolics

###############################################################################
# Semidiscretization of the SWE-Exner equations with source terms for convergence testing

# Equations with Grass model
equations = ShallowWaterExnerEquations1D(gravity = 10.0, rho_f = 0.5,
                                         rho_s = 1.0, porosity = 0.5,
                                         friction = ManningFriction(n = 0.0),
                                         sediment_model = GrassModel(A_g = 0.01))

##################################################################################################
### Create manufactured solution for method of manufactured solutions (MMS)

# Symbolic Variables
@variables x_sym, t_sym, g, r

# Define Differentials
Dt, Dx = Differential(t_sym), Differential(x_sym)

##################################################################################################
##  Initial condition
h = 2 + cos(pi * x_sym) * cos(pi * t_sym)
v = -0.35 * cos(0.5 * pi * x_sym)
h_b = 1 + sin(pi * x_sym) * cos(pi * t_sym)

# directly write in terms of the conservative variables
init = [h, h * v, h_b]

##################################################################################################
##  PDE
# Helper variables for Grass sediment closure
# where we hard code that m_g = 3 to simplify expressions
@variables porosity_inv Ag
h_s = porosity_inv * Ag * v^2
q_s = h_s * v

# PDE Source Terms
eqs = [Dt(h) + Dx(h * v),
    Dt(h * v) + Dx(h * v^2) + g * h * Dx(h + h_b) + g * (h_s / r) * Dx(r * h + h_b),
    Dt(h_b) + Dx(q_s)]

###################################################################################################
## Create the functions for the manufactured solution

function initial_condition_convergence(equations::ShallowWaterExnerEquations1D)
    # Build functions from the symbolic expressions
    init_f1 = build_function(init[1], x_sym, t_sym, expression = Val(false))
    init_f2 = build_function(init[2], x_sym, t_sym, expression = Val(false))
    init_f3 = build_function(init[3], x_sym, t_sym, expression = Val(false))
    function initial_condition_convergence(x, t,
                                           equations::ShallowWaterExnerEquations1D)
        x1 = x[1]
        T = eltype(x)
        return SVector{3, T}(init_f1(x1, t),
                             init_f2(x1, t),
                             init_f3(x1, t))
    end
end

function source_terms_convergence(equations::ShallowWaterExnerEquations1D)
    # Build functions from the symbolic expressions
    du_f1 = build_function(expand_derivatives(eqs[1]),
                           x_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    du_f2 = build_function(expand_derivatives(eqs[2]),
                           x_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    du_f3 = build_function(expand_derivatives(eqs[3]),
                           x_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    return function source_terms_convergence(u, x, t,
                                             equations::ShallowWaterExnerEquations1D)
        g = equations.gravity
        r = equations.r
        porosity_inv = equations.porosity_inv
        Ag = equations.sediment_model.A_g
        x1 = x[1]
        T = eltype(u)
        return SVector{3, T}(du_f1(x1, t, g, r, porosity_inv, Ag),
                             du_f2(x1, t, g, r, porosity_inv, Ag),
                             du_f3(x1, t, g, r, porosity_inv, Ag))
    end
end

initial_condition = initial_condition_convergence(equations)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

solver = DGSEM(polydeg = 4,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = -1.0
coordinates_max = 1.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 2,
                n_cells_max = 10_000,
                periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence(equations),
                                    boundary_conditions = boundary_condition_periodic)

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

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
