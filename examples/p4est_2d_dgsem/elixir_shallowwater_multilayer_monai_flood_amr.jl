
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography

# See the tutorial
# https://trixi-framework.github.io/TrixiShallowWater.jl/stable/tutorials/elixir_shallowwater_monai_tsunami/
# for a thorough explanation of this problem setup.
#
# This elixir allows for experimentation with AMR with this complex test case.
# In the current state, AMR on wet/dry fronts remains sensitive to the limiting of the water height,
# desingularization of the velocities, and the design of the AMR indicator.
# Further, there remains an appreciable loss of conservation, with conservation errors
# of ≈1e-5 in the water height for this test case as the flow evolves.
# More investigation is needed to identify exactly why each aspect of this strange behavior occurs.
# For instance, the conservation loss  may be related to the choice of the CFL that guarantees
# a positive water height as well as a positive element-wise mean necessary in the limiting.

###############################################################################
# Semidiscretization of the multilayer shallow water equations with one layer
# to fallback to be the standard shallow water equations

equations = ShallowWaterMultiLayerEquations2D(gravity = 9.81, H0 = 0.0,
                                              rhos = 1.0)

# Get the precision type for type stable dispatch on the helper spline functions
RealT = typeof(equations.gravity)

# Create a bicubic B-spline interpolation of the bathymetry data, then create a function
# to evaluate the resulting spline at a given point $(x,y)$.
spline_bathymetry_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/21255c980c4eda5294f91e8dfe6c7e33/raw/1afb73928892774dc3a902e0c46ffd882ef03ee3/monai_bathymetry_data.txt",
                                        joinpath(@__DIR__, "monai_bathymetry_data.txt"));

# B-spline interpolation of the underlying data.
# The type of this struct is fixed as `BicubicBSpline`.
const bath_spline_struct = BicubicBSpline(spline_bathymetry_file,
                                          end_condition = "not-a-knot")
bathymetry(x::RealT, y::RealT) = spline_interpolation(bath_spline_struct, x, y)

# Initial condition is a constant background state
function initial_condition_monai_tsunami(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Set the background water height
    H = equations.H0

    # Initially water is at rest
    v1 = zero(H)
    v2 = zero(H)

    # Bottom topography values are computed from the bicubic spline created above
    x1, x2 = x
    b = bathymetry(x1, x2)

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations2D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    H = max(H, b + equations.threshold_limiter)

    # Return the conservative variables
    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_monai_tsunami

# Tsunami test case uses a specialized wave maker type of boundary condition.
# It is used to model an incident wave that approaches from off-shore
# with a water depth of $h = 13.535\,\text{cm}$. To create the incident wave information
# that is valid over the time interval $t \in [0\,s, 22.5\,s]$.

# We download the incident wave data that has been preprocessed to be in the format
# required by TrixiBottomTopography.
water_height_data = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/5b11f5f175bddb326d11d8e28398127e/raw/64980e0e4526e0fcd49589b34ee5458b9a1cebff/monai_wavemaker_bc.txt",
                                   joinpath(@__DIR__, "monai_wavemaker_bc.txt"));

# Similar to the bathymetry approximation, we construct a cubic B-spline interpolation
# of the data, then create a function to evaluate the resulting spline at a given $t$ value.
# The type of this struct is fixed as `CubicBSpline`.
const h_spline_struct = CubicBSpline(water_height_data; end_condition = "not-a-knot")
H_from_wave_maker(t::RealT) = spline_interpolation(h_spline_struct, t)

# Now we can define the specialized boundary condition for the incident wave maker.
@inline function boundary_condition_wave_maker(u_inner, normal_direction::AbstractVector,
                                               x, t, surface_flux_functions,
                                               equations::ShallowWaterMultiLayerEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # Compute the water height from the wave maker input file data
    # and then clip to avoid negative water heights and division by zero
    h_ext = max(equations.threshold_limiter, H_from_wave_maker(t) - u_inner[4])

    # Compute the incoming velocity as in Eq. (10) of the paper
    # - S. Vater, N. Beisiegel, and J. Behrens (2019)
    #   A limiter-based well-balanced discontinuous Galerkin method for shallow-water flows
    #   with wetting and drying: Triangular grids
    #   [DOI: 10.1002/fld.4762](https://doi.org/10.1002/fld.4762)
    h0 = 0.13535 # reference incident water height converted to meters
    v1_ext = 2 * (sqrt(equations.gravity * h_ext) - sqrt(equations.gravity * h0))

    # Create the external solution state in the conservative variables
    u_outer = SVector(h_ext, h_ext * v1_ext, zero(eltype(x)), u_inner[4])

    # Calculate the boundary flux and nonconservative contributions
    flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)

    noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction,
                                                 equations)

    return flux, noncons_flux
end

boundary_condition = Dict(:Bottom => boundary_condition_slip_wall,
                          :Top => boundary_condition_slip_wall,
                          :Right => boundary_condition_slip_wall,
                          :Left => boundary_condition_wave_maker)

# Manning friction source term
@inline function source_terms_manning_friction(u, x, t,
                                               equations::ShallowWaterMultiLayerEquations2D)
    # Grab the conservative variables
    h, hv_1, hv_2, _ = u

    # Desingularization
    h = (h^2 + max(h^2, 1e-8)) / (2 * h)

    # friction coefficient
    n = 0.001

    # Compute the common friction term
    Sf = -equations.gravity * n^2 * h^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    return SVector(zero(h),
                   Sf * hv_1,
                   Sf * hv_2,
                   zero(h))
end

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)

surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(3)

# Cannot simply use `waterheight` here for multilayer equations.
# Need a helper function to extract the relevant variable.
@inline function main_waterheight(u, equations)
    return waterheight(u, equations)[1]
end

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = main_waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the unstructured quad mesh from a file (downloads the file if not available locally)

mesh_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/d26356bdc8e8d742a3035b3f71c71a68/raw/9d6ceedb844e92313d1dac2318a28c87ffbb9de2/mesh_monai_shore.inp",
                           joinpath(@__DIR__, "mesh_monai_shore.inp"));

mesh = P4estMesh{2}(mesh_file)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_terms_manning_friction)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 22.5)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,),
                                     save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.2,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.5)

# Another possible AMR indicator function could be the water height with `IndicatorLoehner`
@inline function momentum_norm(u, equations::ShallowWaterMultiLayerEquations2D)
    _, h_v1, h_v2, _ = u
    return sqrt(h_v1^2 + h_v2^2)
end

amr_controller = ControllerThreeLevel(semi,
                                      IndicatorMax(semi, variable = momentum_norm),
                                      base_level = 0,
                                      med_level = 1, med_threshold = 0.012,
                                      max_level = 3, max_threshold = 0.015)

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 5,
                           adapt_initial_condition = false,
                           adapt_initial_condition_only_refine = false)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        amr_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()...,
            callback = callbacks,
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            adaptive = false);
