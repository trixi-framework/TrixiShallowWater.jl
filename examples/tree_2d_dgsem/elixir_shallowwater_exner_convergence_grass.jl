
OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using Symbolics

###############################################################################
# Semidiscretization of the shallow water Exner equations with source terms for convergence testing

equations = ShallowWaterExnerEquations2D(gravity = 9.81, H0 = 10.0,
                                         rho_f = 0.5, rho_s = 1.0,
                                         porosity = 0.4,
                                         sediment_model = GrassModel(A_g = 0.01))

##################################################################################################
### Create manufactured solution for method of manufactured solutions (MMS)

# Symbolic Variables
@variables x_sym, y_sym, t_sym, g, r

# Define Differentials
Dt, Dx, Dy = Differential(t_sym), Differential(x_sym), Differential(y_sym)

##################################################################################################
##  Initial condition
h = 2 + cos(pi * x_sym) * cos(pi * y_sym) * cos(pi * t_sym)
v1 = 0.5 * cos(0.5 * pi * x_sym) * cos(0.5 * pi * y_sym)
v2 = -0.65 * cos(0.5 * pi * x_sym) * cos(0.5 * pi * y_sym)
h_b = 1 + sin(pi * x_sym) * sin(pi * y_sym) * cos(pi * t_sym)

# directly write in terms of the conservative variables
init = [h, h * v1, h * v2, h_b]

##################################################################################################
##  PDE
# Helper variables for Grass sediment closure
# where we hard code that m_g = 3 to simplify expressions
@variables porosity_inv Ag
h_s = porosity_inv * Ag * (v1^2 + v2^2)
q_s1 = h_s * v1
q_s2 = h_s * v2

# PDE Source Terms
eqs = [Dt(h) + Dx(h * v1) + Dy(h * v2),
    Dt(h * v1) + Dx(h * v1^2) + Dy(h * v1 * v2) + g * h * Dx(h + h_b) +
    g * (h_s / r) * Dx(r * h + h_b),
    Dt(h * v2) + Dx(h * v1 * v2) + Dy(h * v2^2) + g * h * Dy(h + h_b) +
    g * (h_s / r) * Dy(r * h + h_b),
    Dt(h_b) + Dx(q_s1) + Dy(q_s2)]

###################################################################################################
## Create the functions for the manufactured solution

function initial_condition_convergence(equations::ShallowWaterExnerEquations2D)
    # Build functions from the symbolic expressions
    init_f1 = build_function(init[1], x_sym, y_sym, t_sym, expression = Val(false))
    init_f2 = build_function(init[2], x_sym, y_sym, t_sym, expression = Val(false))
    init_f3 = build_function(init[3], x_sym, y_sym, t_sym, expression = Val(false))
    init_f4 = build_function(init[4], x_sym, y_sym, t_sym, expression = Val(false))
    function initial_condition_convergence(x, t,
                                           equations::ShallowWaterExnerEquations2D)
        x1, x2 = x
        T = eltype(x)
        return SVector{4, T}(init_f1(x1, x2, t),
                             init_f2(x1, x2, t),
                             init_f3(x1, x2, t),
                             init_f4(x1, x2, t))
    end
end

function source_terms_convergence(equations::ShallowWaterExnerEquations2D)
    # Build functions from the symbolic expressions
    du_f1 = build_function(expand_derivatives(eqs[1]),
                           x_sym, y_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    du_f2 = build_function(expand_derivatives(eqs[2]),
                           x_sym, y_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    du_f3 = build_function(expand_derivatives(eqs[3]),
                           x_sym, y_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    du_f4 = build_function(expand_derivatives(eqs[4]),
                           x_sym, y_sym, t_sym,
                           g, r, porosity_inv, Ag,
                           expression = Val(false))
    return function source_terms_convergence(u, x, t,
                                             equations::ShallowWaterExnerEquations2D)
        g = equations.gravity
        r = equations.r
        porosity_inv = equations.porosity_inv
        Ag = equations.sediment_model.A_g
        x1, x2 = x
        T = eltype(u)
        return SVector{4, T}(du_f1(x1, x2, t, g, r, porosity_inv, Ag),
                             du_f2(x1, x2, t, g, r, porosity_inv, Ag),
                             du_f3(x1, x2, t, g, r, porosity_inv, Ag),
                             du_f4(x1, x2, t, g, r, porosity_inv, Ag))
    end
end

initial_condition = initial_condition_convergence(equations)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxPlusDissipation(flux_ersing_etal, dissipation_roe),
                flux_nonconservative_ersing_etal)
solver = DGSEM(polydeg = 3,
               surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_convergence(equations),
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver

tspan = (0.0, 0.75)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback)

###############################################################################
# run the simulation
sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
