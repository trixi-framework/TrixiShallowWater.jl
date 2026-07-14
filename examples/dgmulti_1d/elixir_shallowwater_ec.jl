using OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the shallow water equations with a discontinuous
# bottom topography function for a fully wet configuration

equations = ShallowWaterEquations1D(gravity = 9.81)

# Initial condition with a truly discontinuous water height, velocity, and bottom
# topography function as an academic testcase for entropy conservation.
# The errors from the analysis callback are not important but `∑∂S/∂U ⋅ Uₜ` should
# be around machine roundoff.
# Works as intended for a `DGMultiMesh` with `cells_per_dimension = (16,)`. If the mesh
# resolution is changed the initial condition below may need changed as well to
# ensure that the discontinuities lie on an element interface.
function initial_condition_ec_discontinuous_bottom(x, t,
                                                   equations::ShallowWaterEquations1D)
    # Set the background values
    H = 4.25
    v = 0.0
    b = sin(x[1]) # arbitrary continuous function

    # Setup the discontinuous water height and velocity
    if x[1] >= 0.125 && x[1] <= 0.25
        H = 5.0
        v = 0.1882
    end

    # Setup a discontinuous bottom topography
    if x[1] >= -0.25 && x[1] <= -0.125
        b = 2.0 + 0.5 * sin(2.0 * pi * x[1])
    end

    return prim2cons(SVector(H, v, b), equations)
end

initial_condition = initial_condition_ec_discontinuous_bottom

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal)

dg = DGMulti(polydeg = 4, element_type = Line(), approximation_type = SBP(),
             surface_integral = SurfaceIntegralWeakForm(surface_flux),
             volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the DGMulti mesh and setup a periodic mesh

cells_per_dimension = (16,)
mesh = DGMultiMesh(dg, cells_per_dimension,
                   coordinates_min = (-1.0,), coordinates_max = (1.0,),
                   periodicity = true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg;
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# ODE solver

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     uEltype = real(dg))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation

sol = solve(ode, RDPK3SpFSAL49(); abstol = 1.0e-8, reltol = 1.0e-8,
            ode_default_options()..., callback = callbacks)
