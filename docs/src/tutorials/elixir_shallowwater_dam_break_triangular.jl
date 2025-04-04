# Note: This tutorial is still under construction

# In this tutorial we will use the shallow water equations to simulate a dam break over triangular 
# bottom topography with wetting and drying.
# More information about this test can be found in the papers:
# - S. Gu et al. (2017)
#   SWE-SPHysics Simulation of Dam Break Flows at South-Gate Gorges Reservoir
#   [DOI: 10.3390/w9060387](https://doi.org/10.3390/w9060387)
# - J.G. Zhou, D.M. Causon, C.G. Mingham and D.M. Ingram (2004)
#   Numerical Prediction of Dam-Break Flows in General Geometries with Complex Bed Topography
#   [DOI: 10.1061/(ASCE)0733-9429(2004)130:4(332)](https://doi.org/10.1061/(ASCE)0733-9429(2004)130:4(332))

# The tutorial will cover:
# - Set up a SWE solver for wet/dry transitions
# - Create custom initial conditions / source terms
# - save solution data at gauge points
# - Visualization

# Before we start, we need to load the required packages. Besides TrixiShallowWater.jl, we require
# [`Trixi.jl`](@extref Trixi.jl) for the spatial discretization and `OrdinaryDiffEqSSPRK.jl` for time integration,
# and `CairoMakie.jl` for visualization.
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using CairoMakie

# In the first step we will set up the equation system. In this example we want to solve the 
# one-dimensional shallow water equations, so we will use the [`ShallowWaterEquationsWetDry1D`](@ref ShallowWaterEquationsWetDry1D)
# and specify the gravitational acceleration to `gravity_constant = 9.81`. In contrast to the
# [`Trixi.ShallowWaterEquations1D`](@extref Trixi.ShallowWaterEquations1D) type, this equation type 
# allows contains additional parameters and methods needed to handle wetting and drying.
equations = ShallowWaterEquationsWetDry1D(gravity_constant = 9.812)

# We then create a function to supply the initial condition for the simulation.
function initial_condition_dam_break_triangular(x, t,
                                                equations::ShallowWaterEquationsWetDry1D)
    b = 0.0  # Bottom topography
    h = 0.0  # Water height
    v = 0.0  # Velocity

    if x[1] <= 15.5
        h = 0.75  # Water height in the left reservoir
    elseif 25.5 < x[1] && x[1] <= 28.5
        b = (x[1] - 25.5) * 0.4 / 3.0  # Rising slope of the triangular bottom
    elseif x[1] > 28.5 && x[1] < 31.5
        b = 0.4 - (x[1] - 28.5) * 0.4 / 3.0  # Falling slope of the triangular bottom
    end

    H = h + b  # Total water height
    if x[1] > 28.5
        H = max(H, 0.15)  # Water height in the right reservor
    end

    ## It is mandatory to shift the water level at dry areas to make sure the water height h
    ## stays positive. The system would not be stable for h set to a hard zero due to division by h in
    ## the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    ## with a default value of 5*eps() ≈ 1e-13 in double precision, is set in the constructor above
    ## for the ShallowWaterEquations and added to the initial condition if h = 0.
    ## This default value can be changed within the constructor call depending on the simulation setup.
    H = max(H, b + equations.threshold_limiter)
    return prim2cons(SVector(H, v, b), equations)
end

initial_condition = initial_condition_dam_break_triangular;

# Let's now use `Makie.jl` to visualize the setup.
x_init = collect(0:0.1:38.0)
fig, ax, plt = lines(x_init,
                     [cons2prim(initial_condition(x, 0.0, equations), equations)[1]
                      for x in x_init], color = :blue)
lines!(ax, x_init, [initial_condition(x, 0.0, equations)[3] for x in x_init],
       color = :black, linestyle = :solid)
fig

# As we want to compare the results to experimental data, we also need to model the bottom friction.
# For this we create a new source term, which adds a Manning friction term to the momentum equations.
@inline function source_term_manning_friction(u, x, t,
                                              equations::ShallowWaterEquationsWetDry1D)
    h, hv, _ = u

    n = 0.0125  # friction coefficient
    h = (h^2 + max(h^2, 1e-8)) / (2.0 * h) # desingularization procedure

    return SVector(0.0, -equations.gravity * n^2 * h^(-7 / 3) * abs(hv) * hv, 0.0)
end

# Now we can set up the approximation space, where we use the discontinuous Galerkin spectral element
# method ([`DGSEM`](@extref Trixi.DGSEM)), with a volume integral in flux differencing formulation. 
# For this we first need to specify fluxes for both volume and surface integrals. Since the system 
# is setup in nonconservative form the fluxes need to provided in form of a tuple 
# `flux = (conservative flux, nonconservative_flux)`. To ensure well-balancedness and positivity a 
# reconstruction procedure is applied for the surface fluxes and a special shock-capturing scheme
# is used to compute the volume integrals.
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(4)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = Trixi.waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

# The mesh is created using the `TreeMesh` type. 
# The computational domain spans from `coordinates_min` to `coordinates_max` and is initialized with
# 2^7 = 128 elements. The domain is non-periodic.
coordinates_min = 0.0
coordinates_max = 38.0

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 8,
                n_cells_max = 10_000,
                periodicity = false)

# The semi-discretization object combines the mesh, equations, initial condition,
# solver, boundary conditions, and source terms into a single object. This object
# represents the spatial discretization of the problem and is complemented with the required 
# time interval to define an ODE problem for time integration.
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition_slip_wall,
                                    source_terms = source_term_manning_friction)
tspan = (0.0, 40.0)
ode = semidiscretize(semi, tspan);

# Callbacks are used to monitor the simulation, save results, and control the time step size.
# Below, we define several callbacks for different purposes.

# ### Analysis Callback
# Performs analysis at regular intervals, such as computing conservation errors.
analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,),
                                     save_analysis = false)

# ### Time Series Callback
# Saves time series data at specific gauge locations for further analysis.
time_series = TimeSeriesCallback(semi, [(19.5), (25.5), (28.5), (35.5)];
                                 interval = 1,
                                 solution_variables = cons2cons,
                                 filename = "tseries.h5")

# ### Stepsize Callback
# Controls the time step size based on the CFL condition.
stepsize_callback = StepsizeCallback(cfl = 0.5)

# ### Combine Callbacks
# All the defined callbacks are combined into a single `CallbackSet`.
callbacks = CallbackSet(analysis_callback,
                        time_series,
                        stepsize_callback);

# ## Running the Simulation
#
# Finally, we solve the ODE problem using a strong stability-preserving Runge-Kutta (SSPRK) method.
# The `PositivityPreservingLimiterShallowWater` is used as a stage limiter to ensure positivity
# of the water height during the simulation. The `SSPRK43` integrator supports adaptive timestepping,
# but since we use a CFL-based time step we set (`adaptive = false`).
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))
save_times = (0, 3.5, 40)
sol = solve(ode, SSPRK43(stage_limiter!); dt = 1.0,
            ode_default_options()..., callback = callbacks, adaptive = false,
            saveat = save_times);

# ## Visualization
pd_list = [PlotData1D(sol.u[i], semi, reinterpolate = false) for i in 1:length(save_times)]

f = Figure(size = (550, 550 / 2.5))
ax = Axis(f[1, 1], xlabel = "x / m", ylabel = "waterheight / m", limits = (0, 38, 0.0, 1.2))
for (i, pd) in enumerate(pd_list)
    lines!(ax, pd.x, pd.data[:, 1], label = "t = $(sol.t[i])s")
end
lines!(ax, pd_list[1].x, pd_list[1].data[:, 3], color = :black, linestyle = :solid)
band!(ax, pd_list[1].x, 0.0, pd_list[1].data[:, 3], color = :gray95)  # Set color for bottom topography
axislegend(ax, orientation = :horizontal)
f

# TODO: Download reference data
# Time series data
pd = PlotData1D(time_series, 1)
lines(pd.x, pd.data[:, 1])
