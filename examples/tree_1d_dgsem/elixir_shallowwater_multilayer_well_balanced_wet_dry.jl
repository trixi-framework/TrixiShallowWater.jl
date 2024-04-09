
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the two-layer shallow water equations to test well-balancedness

equations = ShallowWaterMultiLayerEquations1D(gravity_constant = 1.0, H0 = 2.0, rhos = (1.0, 3.0))

"""
    initial_condition_twolayer_well_balanced_wet_dry(perturbation, equations:: ShallowWaterMultiLayerEquations1D)
    initial_condition_twolayer_well_balanced_wet_dry(x, t, equations:: ShallowWaterMultiLayerEquations1D)

Initial condition with a complex (discontinuous) bottom topography to test well-balancedness for a
two-layer shallow water system with dry states if `perturbation` is set to `false`.
Additionally, it is possible to set a perturbation in the lower layer to test perturbations from the 
lake-at-rest condition by setting the `perturbation` variable to `true`.

The initial condition is taken from the paper:
  - S. Martínez-Aranda, A. Ramos-Pérez, P. García-Navarro (2020)
    A 1D shallow-ﬂow model for two-layer ﬂows based on FORCE scheme with wet–dry treatment\n 
    [DOI:10.2166/hydro.2020.002](https://doi.org/10.2166/hydro.2020.002)
"""
function initial_condition_twolayer_well_balanced_wet_dry(perturbation::Bool, equations::ShallowWaterMultiLayerEquations1D)
    return function initial_condition_twolayer_well_balanced_wet_dry(x, t, equations)
        # Set interface height for the upper layer
        H_upper = 2.0

        # Set interface height for the lower layer
        H_lower = 1.0
        if perturbation == true
            if x[1] <= 5.0
                nothing
            elseif x[1] <= 15.0
                H_lower = 0.75
            elseif x[1] <= 35.0
                H_lower = 1.5 - (x[1] - 15.0) / 20.0 * 0.5
            elseif x[1] <= 40.0
                H_lower = 1.0 + (x[1] - 35.0) / 5.0 * 0.5
            elseif x[1] <= 55.0
                H_lower = 1.0
            elseif x[1] <= 70.0
                H_lower = 1.6 - (x[1] - 55.0) / 15.0 * 0.4
            elseif x[1] <= 75.0
                nothing
            elseif x[1] <= 85.0
                H_lower = 1.0 - (x[1] - 75.0) / 10.0 * 0.3
            elseif x[1] <= 95.0
                H_lower = 1.1
            else
                H_lower = 1.5
            end
        end

        # Set the bottom topography
        if x[1] <= 5.0
            b = 2.5
        elseif x[1] <= 20.0
            b = 0.5
        elseif x[1] <= 40.0
            b = 0.1 + (x[1] - 20.0) / 20.0 * 1.4
        elseif x[1] <= 60.0
            b = 1.0 - (x[1] - 40.0) / 20.0 * 1.3
        elseif x[1] <= 70.0
            b = 0.7
        elseif x[1] <= 75.0
            b = 2.2
        elseif x[1] <= 95.0
            b = 0.4
        else
            b = 1.5
        end

        # Set zero initial velocity
        v1_upper = 0.0
        v1_lower = 0.0

        #= 
        It is mandatory to shift the water level at dry areas to make sure the water height h
        stays positive. The system would not be stable for h set to a hard 0 due to division by h in
        the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
        with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
        for the ShallowWaterMultiLayerEquations1D and added to the initial condition if h = 0.
        This default value can be changed within the constructor call depending on the simulation setup.
        =#
        H_lower = max(H_lower, b + equations.threshold_limiter)
        H_upper = max(H_upper, H_lower + equations.threshold_limiter)

        return prim2cons(SVector(H_upper, H_lower, v1_upper, v1_lower, b), equations)
    end
end


initial_condition = initial_condition_twolayer_well_balanced_wet_dry(true, equations)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))
basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)
###############################################################################
# Get the TreeMesh and setup a periodic mesh

# The domain of interest [0.0, 100.0] is extended towards -28.0 to match the positions of the 
# disctontinuities with TreeMesh.
coordinates_min = -28.0
coordinates_max = 100.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 10_000,
                periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

#################################################################################
# Workaround to set a discontinuous water and bottom topography for
# debugging and testing. Essentially, this is a slight augmentation of the
# `compute_coefficients` where the `x` node value passed here is slightly
# perturbed to the left / right in order to set a true discontinuity that avoids
# the doubled value of the LGL nodes at a particular element interface.
#
# Note! The errors from the analysis callback are not important but the error
# for this lake at rest test case `∑|H0-(h+b)|` should be near machine roundoff.

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for i in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, element)
        # We know that the discontinuity is a vertical line. Slightly augment the x value by a factor
        # of unit roundoff to avoid the repeted value from the LGL nodes at at interface.
        if i == 1
            x_node = SVector(nextfloat(x_node[1]))
        elseif i == nnodes(semi.solver)
            x_node = SVector(prevfloat(x_node[1]))
        end
        u_node = initial_condition(x_node, first(tspan), equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, element)
    end
end
#############################################################################################

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_integrals = (lake_at_rest_error,))

stepsize_callback = StepsizeCallback(cfl = 0.5)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

# use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);
summary_callback() # print the timer summary
