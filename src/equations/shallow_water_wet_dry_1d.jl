# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    ShallowWaterEquationsWetDry1D(; gravity, H0 = 0, threshold_limiter = nothing threshold_wet = nothing)

Shallow water equations (SWE) in one space dimension. The equations are given by
```math
\begin{aligned}
  \frac{\partial h}{\partial t} + \frac{\partial}{\partial x}(h v) &= 0 \\
    \frac{\partial}{\partial t}(h v) + \frac{\partial}{\partial x}\left(h v^2 + \frac{g}{2}h^2\right)
    + g h \frac{\partial b}{\partial x} &= 0
\end{aligned}
```
The unknown quantities of the SWE are the water height ``h`` and the velocity ``v``.
The gravitational acceleration is denoted by `g` and the (possibly) variable bottom topography
function ``b(x)``. Conservative variable water height ``h`` is measured from the bottom topography ``b``, therefore one
also defines the total water height as ``H = h + b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water height that
is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

Also, there are four thresholds which prevent numerical problems as well as instabilities. None of them
have to be passed, as default values are defined within the struct. The first one, `threshold_limiter`, is
used in [`PositivityPreservingLimiterShallowWater`](@ref) on the water height, as a (small) shift on the initial
condition and cutoff before the next time step. The second one, `threshold_wet`, is applied on the water height to
define when the flow is "wet" before calculating the numerical flux. A third
`threshold_partially_wet` is applied on the water height to define "partially wet" elements in
[`IndicatorHennemannGassnerShallowWater`](@ref), that are then calculated with a pure FV method to
ensure well-balancedness. Lastly, `threshold_desingularization` is used in [`PositivityPreservingLimiterShallowWater`](@ref)
for the velocity desingularization procedure. For `Float64` no threshold needs to be passed, as default values are
defined within the struct. For other number formats `threshold_partially_wet`
and `threshold_desingularization` must be provided.

The bottom topography function ``b(x)`` is set inside the initial condition routine
for a particular problem setup. To test the conservative form of the SWE one can set the bottom topography
variable `b` to zero.

In addition to the unknowns, Trixi.jl currently stores the bottom topography values at the approximation points
despite being fixed in time. This is done for convenience of computing the bottom topography gradients
on the fly during the approximation as well as computing auxiliary quantities like the total water height ``H``
or the entropy variables.
This affects the implementation and use of these equations in various ways:
* The flux values corresponding to the bottom topography must be zero.
* The bottom topography values must be included when defining initial conditions, boundary conditions or
  source terms.
* [`Trixi.AnalysisCallback`](@extref) analyzes this variable.
* Trixi.jl's visualization tools will visualize the bottom topography by default.

References for the SWE are many but a good introduction is available in Chapter 13 of the book:
- Randall J. LeVeque (2002)
  Finite Volume Methods for Hyperbolic Problems
  [DOI: 10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253)
"""
struct ShallowWaterEquationsWetDry1D{RealT <: Real} <:
       Trixi.AbstractShallowWaterEquations{1, 3}
    gravity::RealT # gravitational acceleration
    H0::RealT      # constant "lake-at-rest" total water height
    # `threshold_limiter` used in `PositivityPreservingLimiterShallowWater` on water height,
    # as a (small) shift on the initial condition and cutoff before the next time step.
    # Default is 500*eps() which in double precision is ≈1e-13.
    threshold_limiter::RealT
    # `threshold_wet` applied on water height to define when the flow is "wet"
    # before calculating the numerical flux.
    # Default is 5*eps() which in double precision is ≈1e-15.
    threshold_wet::RealT
    # `threshold_partially_wet` used in `IndicatorHennemannGassnerShallowWater` on the water height
    # to define "partially wet" elements. Those elements are calculated with a pure FV method to
    # ensure well-balancedness. Default in double precision is 1e-4.
    threshold_partially_wet::RealT
    # `threshold_desingularization` used in the velocity desingularization procedure, to avoid
    # division by small numbers. Default in double precision is 1e-10.
    threshold_desingularization::RealT
    # Standard shallow water equations for dispatch on Trixi.jl functions
    basic_swe::ShallowWaterEquations1D{RealT}
end

# Allow for flexibility to set the gravitational acceleration within an elixir depending on the
# application where `gravity=1.0` or `gravity=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases.
# Strict default values for thresholds that performed well in many numerical experiments
function ShallowWaterEquationsWetDry1D(; gravity, H0 = zero(gravity),
                                       threshold_limiter = nothing,
                                       threshold_wet = nothing,
                                       threshold_partially_wet = nothing,
                                       threshold_desingularization = nothing)
    T = promote_type(typeof(gravity), typeof(H0))
    if threshold_limiter === nothing
        threshold_limiter = 500 * eps(T)
    end
    if threshold_wet === nothing
        threshold_wet = 5 * eps(T)
    end
    if threshold_partially_wet === nothing
        threshold_partially_wet = default_threshold_partially_wet(T)
    end
    if threshold_desingularization === nothing
        threshold_desingularization = default_threshold_desingularization(T)
    end
    # Construct the standard SWE for dispatch. Even though the `basic_swe` already store the
    # gravitational acceleration and the total water height, we store an extra copy in
    # `ShallowWaterEquationsWetDry1D` for convenience.
    basic_swe = ShallowWaterEquations1D(gravity, H0)

    ShallowWaterEquationsWetDry1D(gravity, H0, threshold_limiter,
                                  threshold_wet, threshold_partially_wet,
                                  threshold_desingularization, basic_swe)
end

Trixi.have_nonconservative_terms(::ShallowWaterEquationsWetDry1D) = True()

function Trixi.varnames(::typeof(cons2cons), ::ShallowWaterEquationsWetDry1D)
    ("h", "h_v", "b")
end
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
function Trixi.varnames(::typeof(cons2prim), ::ShallowWaterEquationsWetDry1D)
    ("H", "v", "b")
end

# This equation set extends the basic ShallowWaterEquations1D from Trixi.jl with additional functionality
# for wet/dry transitions. Since many functions correspond to the fully wet case, we make use of the
# existing functionality and introduce a number of wrapper functions, that dispatch to the
# ShallowWaterEquations1D.

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterEquationsWetDry1D)

A smooth initial condition used for convergence tests in combination with
[`source_terms_convergence_test`](@ref)
(and [`Trixi.BoundaryConditionDirichlet`](@extref) in non-periodic domains).
"""
function Trixi.initial_condition_convergence_test(x, t,
                                                  equations::ShallowWaterEquationsWetDry1D)
    return Trixi.initial_condition_convergence_test(x, t,
                                                    equations.basic_swe)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterEquationsWetDry1D)

Source terms used for convergence tests in combination with
[`initial_condition_convergence_test`](@ref)
(and [`Trixi.BoundaryConditionDirichlet`](@extref) in non-periodic domains).

This manufactured solution source term is specifically designed for the bottom topography function
`b(x) = 2.0 + 0.5 * sin(sqrt(2.0)*pi*x[1])`
as defined in [`initial_condition_convergence_test`](@ref).
"""

@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterEquationsWetDry1D)
    return Trixi.source_terms_convergence_test(u, x, t,
                                               equations.basic_swe)
end

"""
    initial_condition_weak_blast_wave(x, t, equations::ShallowWaterEquationsWetDry1D)

A weak blast wave discontinuity useful for testing, e.g., total energy conservation.
Note for the shallow water equations to the total energy acts as a mathematical entropy function.
"""
function Trixi.initial_condition_weak_blast_wave(x, t,
                                                 equations::ShallowWaterEquationsWetDry1D)
    return Trixi.initial_condition_weak_blast_wave(x, t,
                                                   equations.basic_swe)
end

"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, x, t, surface_flux_function,
                                  equations::ShallowWaterEquationsWetDry1D)

Create a boundary state by reflecting the normal velocity component and keep
the tangential velocity component unchanged. The boundary water height is taken from
the internal value.

For details see Section 9.2.5 of the book:
- Eleuterio F. Toro (2001)
  Shock-Capturing Methods for Free-Surface Shallow Flows
  1st edition
  ISBN 0471987662
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation_or_normal,
                                                    direction,
                                                    x, t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterEquationsWetDry1D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # This can not be dispatched, when Flux Hydrostactic reconstruction is used
    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         -u_inner[2],
                         u_inner[3])

    # calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal,
                                     equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary,
                                                     orientation_or_normal,
                                                     equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal,
                                     equations)
        noncons_flux = nonconservative_flux_function(u_boundary, u_inner,
                                                     orientation_or_normal,
                                                     equations)
    end

    return flux, noncons_flux
end

# Calculate 1D flux for a single point
# Note, the bottom topography has no flux
@inline function Trixi.flux(u, orientation::Integer,
                            equations::ShallowWaterEquationsWetDry1D)
    return Trixi.flux(u, orientation,
                      equations.basic_swe)
end

"""
    flux_nonconservative_wintermeyer_etal(u_ll, u_rr, orientation::Integer,
                                          equations::ShallowWaterEquationsWetDry1D)

Non-symmetric two-point volume flux discretizing the nonconservative (source) term
that contains the gradient of the bottom topography [`ShallowWaterEquationsWetDry1D`](@ref).

Gives entropy conservation and well-balancedness on both the volume and surface when combined with
[`flux_wintermeyer_etal`](@ref).

Further details are available in the papers:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on
  curvilinear meshes
  [DOI: 10.48550/arXiv.2306.12699](https://doi.org/10.48550/arXiv.2306.12699)
"""
@inline function Trixi.flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                             orientation::Integer,
                                                             equations::ShallowWaterEquationsWetDry1D)
    return Trixi.flux_nonconservative_wintermeyer_etal(u_ll, u_rr, orientation,
                                                       equations.basic_swe)
end

"""
    flux_nonconservative_fjordholm_etal(u_ll, u_rr, orientation::Integer,
                                        equations::ShallowWaterEquationsWetDry1D)

Non-symmetric two-point surface flux discretizing the nonconservative (source) term of
that contains the gradient of the bottom topography [`ShallowWaterEquationsWetDry1D`](@ref).

This flux can be used together with [`flux_fjordholm_etal`](@ref) at interfaces to ensure entropy
conservation and well-balancedness.

Further details for the original finite volume formulation are available in
- Ulrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011)
  Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography
  [DOI: 10.1016/j.jcp.2011.03.042](https://doi.org/10.1016/j.jcp.2011.03.042)
and for curvilinear 2D case in the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
"""
@inline function Trixi.flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                           orientation::Integer,
                                                           equations::ShallowWaterEquationsWetDry1D)
    return Trixi.flux_nonconservative_fjordholm_etal(u_ll, u_rr, orientation,
                                                     equations.basic_swe)
end

"""
    flux_nonconservative_audusse_etal(u_ll, u_rr, orientation::Integer,
                                      equations::ShallowWaterEquationsWetDry1D)

Non-symmetric two-point surface flux that discretizes the nonconservative (source) term.
The discretization uses the [`hydrostatic_reconstruction_audusse_etal`](@ref) on the conservative
variables.

This hydrostatic reconstruction ensures that the finite volume numerical fluxes remain
well-balanced for discontinuous bottom topographies [`ShallowWaterEquationsWetDry1D`](@ref).
Should be used together with [`Trixi.FluxHydrostaticReconstruction`](@extref) and
[`hydrostatic_reconstruction_audusse_etal`](@ref) in the surface flux to ensure consistency.

Further details on the hydrostatic reconstruction and its motivation can be found in
- Emmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004)
  A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
  [DOI: 10.1137/S1064827503431090](https://doi.org/10.1137/S1064827503431090)
"""
@inline function Trixi.flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                         orientation::Integer,
                                                         equations::ShallowWaterEquationsWetDry1D)
    return Trixi.flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                   orientation,
                                                   equations.basic_swe)
end

"""
    flux_nonconservative_chen_noelle(u_ll, u_rr,
                                     orientation::Integer,
                                     equations::ShallowWaterEquationsWetDry1D)

Non-symmetric two-point surface flux that discretizes the nonconservative (source) term.
The discretization uses the [`hydrostatic_reconstruction_chen_noelle`](@ref) on the conservative
variables.

Should be used together with [`Trixi.FluxHydrostaticReconstruction`](@extref) and
[`hydrostatic_reconstruction_chen_noelle`](@ref) in the surface flux to ensure consistency.

Further details on the hydrostatic reconstruction and its motivation can be found in
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI:10.1137/15M1053074](https://dx.doi.org/10.1137/15M1053074)
"""
@inline function flux_nonconservative_chen_noelle(u_ll, u_rr,
                                                  orientation::Integer,
                                                  equations::ShallowWaterEquationsWetDry1D)

    # Pull the water height and bottom topography on the left
    h_ll, _, b_ll = u_ll
    h_rr, _, b_rr = u_rr

    H_ll = h_ll + b_ll
    H_rr = h_rr + b_rr

    b_star = min(max(b_ll, b_rr), min(H_ll, H_rr))

    # Create the hydrostatic reconstruction for the left solution state
    u_ll_star, _ = hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, equations)

    # Copy the reconstructed water height for easier to read code
    h_ll_star = u_ll_star[1]

    z = zero(eltype(u_ll))

    return SVector(z,
                   -equations.gravity * (h_ll_star + h_ll) * (b_ll - b_star),
                   z)
end

"""
    flux_fjordholm_etal(u_ll, u_rr, orientation,
                        equations::ShallowWaterEquationsWetDry1D)

Total energy conservative (mathematical entropy for shallow water equations). When the bottom topography
is nonzero this should only be used as a surface flux otherwise the scheme will not be well-balanced.
For well-balancedness in the volume flux use [`flux_wintermeyer_etal`](@ref).

Details are available in Eq. (4.1) in the paper:
- Ulrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011)
  Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography
  [DOI: 10.1016/j.jcp.2011.03.042](https://doi.org/10.1016/j.jcp.2011.03.042)
"""
@inline function Trixi.flux_fjordholm_etal(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquationsWetDry1D)
    return Trixi.flux_fjordholm_etal(u_ll, u_rr, orientation,
                                     equations.basic_swe)
end

"""
    flux_wintermeyer_etal(u_ll, u_rr, orientation,
                          equations::ShallowWaterEquationsWetDry1D)

Total energy conservative (mathematical entropy for shallow water equations) split form.
When the bottom topography is nonzero this scheme will be well-balanced when used as a `volume_flux`.
For the `surface_flux` either [`flux_wintermeyer_etal`](@ref) or [`flux_fjordholm_etal`](@ref) can
be used to ensure well-balancedness and entropy conservation.

Further details are available in Theorem 1 of the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
"""
@inline function Trixi.flux_wintermeyer_etal(u_ll, u_rr, orientation::Integer,
                                             equations::ShallowWaterEquationsWetDry1D)
    return Trixi.flux_wintermeyer_etal(u_ll, u_rr, orientation,
                                       equations.basic_swe)
end

"""
    hydrostatic_reconstruction_audusse_etal(u_ll, u_rr, orientation::Integer,
                                            equations::ShallowWaterEquationsWetDry1D)

A particular type of hydrostatic reconstruction on the water height to guarantee well-balancedness
for a general bottom topography [`ShallowWaterEquationsWetDry1D`](@ref). The reconstructed solution states
`u_ll_star` and `u_rr_star` variables are then used to evaluate the surface numerical flux at the interface.
Use in combination with the generic numerical flux routine [`Trixi.FluxHydrostaticReconstruction`](@extref).

Further details on this hydrostatic reconstruction and its motivation can be found in
- Emmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004)
  A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
  [DOI: 10.1137/S1064827503431090](https://doi.org/10.1137/S1064827503431090)
"""
@inline function Trixi.hydrostatic_reconstruction_audusse_etal(u_ll, u_rr,
                                                               equations::ShallowWaterEquationsWetDry1D)
    return Trixi.hydrostatic_reconstruction_audusse_etal(u_ll, u_rr,
                                                         equations.basic_swe)
end

"""
    hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquationsWetDry1D)

A particular type of hydrostatic reconstruction of the water height to guarantee well-balancedness
for a general bottom topography of the [`ShallowWaterEquationsWetDry1D`](@ref). The reconstructed solution states
`u_ll_star` and `u_rr_star` variables are used to evaluate the surface numerical flux at the interface.
The key idea is a linear reconstruction of the bottom and water height at the interfaces using subcells.
Use in combination with the generic numerical flux routine [`Trixi.FluxHydrostaticReconstruction`](@extref).

Further details on this hydrostatic reconstruction and its motivation can be found in
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI:10.1137/15M1053074](https://dx.doi.org/10.1137/15M1053074)
"""
@inline function hydrostatic_reconstruction_chen_noelle(u_ll, u_rr,
                                                        equations::ShallowWaterEquationsWetDry1D)
    # Unpack left and right water heights and bottom topographies
    h_ll, _, b_ll = u_ll
    h_rr, _, b_rr = u_rr

    # Get the velocities on either side
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    H_ll = b_ll + h_ll
    H_rr = b_rr + h_rr

    b_star = min(max(b_ll, b_rr), min(H_ll, H_rr))

    # Compute the reconstructed water heights
    h_ll_star = min(H_ll - b_star, h_ll)
    h_rr_star = min(H_rr - b_star, h_rr)

    # Set the water height to be at least the value stored in the variable threshold after
    # the hydrostatic reconstruction is applied and before the numerical flux is calculated
    # to avoid numerical problem with arbitrary small values. Interfaces with a water height
    # lower or equal to the threshold can be declared as dry.
    # The default value for `threshold_wet` is ≈ 5*eps(), or 1e-15 in double precision, is set
    # in the `ShallowWaterEquationsWetDry1D` struct. This threshold value can be changed in the constructor
    # call of this equation struct in an elixir.
    threshold = equations.threshold_wet

    if (h_ll_star <= threshold)
        h_ll_star = threshold
        v_ll = zero(v_ll)
    end

    if (h_rr_star <= threshold)
        h_rr_star = threshold
        v_rr = zero(v_rr)
    end

    # Create the conservative variables using the reconstruted water heights
    u_ll_star = SVector(h_ll_star, h_ll_star * v_ll, b_ll)
    u_rr_star = SVector(h_rr_star, h_rr_star * v_rr, b_rr)

    return u_ll_star, u_rr_star
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the bottom topography
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::ShallowWaterEquationsWetDry1D)
    return (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                        orientation_or_normal_direction,
                                                        equations.basic_swe)
end

# Specialized `FluxHLL` to avoid spurious dissipation in the bottom topography
@inline function (numflux::FluxHLL)(u_ll, u_rr, orientation_or_normal_direction,
                                    equations::ShallowWaterEquationsWetDry1D)
    λ_min, λ_max = numflux.min_max_speed(u_ll, u_rr, orientation_or_normal_direction,
                                         equations)

    if λ_min >= 0 && λ_max >= 0
        return flux(u_ll, orientation_or_normal_direction, equations)
    elseif λ_max <= 0 && λ_min <= 0
        return flux(u_rr, orientation_or_normal_direction, equations)
    else
        f_ll = flux(u_ll, orientation_or_normal_direction, equations)
        f_rr = flux(u_rr, orientation_or_normal_direction, equations)
        inv_λ_max_minus_λ_min = inv(λ_max - λ_min)
        factor_ll = λ_max * inv_λ_max_minus_λ_min
        factor_rr = λ_min * inv_λ_max_minus_λ_min
        factor_diss = λ_min * λ_max * inv_λ_max_minus_λ_min
        diss = u_rr - u_ll
        return factor_ll * f_ll - factor_rr * f_rr +
               factor_diss * SVector(diss[1], diss[2], zero(eltype(u_ll)))
    end
end

"""
    dissipation_roe(u_ll, u_rr, orientation_or_normal_direction,
                                    equations::ShallowWaterEquationsWetDry1D)
Roe-type dissipation term for the [`ShallowWaterEquationsWetDry1D`](@ref). To create the classical Roe solver,
this dissipation term can be combined with [`Trixi.flux_central`](@extref) using [`Trixi.FluxPlusDissipation`](@extref).

For details on the Roe linearization see Chapter 15.3.2 and Chapter 15.3.3 for the one-dimensional
shallow water equations of the book:
- Randall J. LeVeque (2002)
  Finite Volume Methods for Hyperbolic Problems
  [DOI: 10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253)
"""
@inline function dissipation_roe(u_ll, u_rr, orientation_or_normal_direction,
                                 equations::ShallowWaterEquationsWetDry1D)
    g = equations.gravity
    z = zero(eltype(u_ll))

    # Get velocities and waterheights
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Compute Roe averages
    h_avg = 0.5f0 * (h_ll + h_rr)
    v_avg = (sqrt(h_ll) * v_ll + sqrt(h_rr) * v_rr) /
            (sqrt(h_ll) + sqrt(h_rr))
    c_avg = (sqrt(g * h_avg))

    # Compute the eigenvalues
    λ1 = v_avg - c_avg
    λ2 = v_avg + c_avg

    # Eigenvector matrix
    R = @SMatrix [[1 1]; [(v_avg - c_avg) (v_avg + c_avg)]]

    # Inverse eigenvector matrix
    R_inv = 1 / (2 * c_avg) * @SMatrix [[(v_avg + c_avg) -1]; [(c_avg - v_avg) 1]]

    # Eigenvalue absolute value matrix
    Λ_abs = @SMatrix [[abs(λ1) z]; [z abs(λ2)]]

    # Compute the jump in conserved variables, excluding the bottom topography
    u_jump = @views (u_rr - u_ll)[1:2]

    diss = SVector(-0.5f0 * R * Λ_abs * R_inv * u_jump)

    return SVector(diss[1], diss[2], z)
end

"""
    min_max_speed_chen_noelle(u_ll, u_rr, orientation::Integer,
                              equations::ShallowWaterEquations1D)

The approximated speeds for the HLL type numerical flux used by Chen and Noelle for their
hydrostatic reconstruction. As they state in the paper, these speeds are chosen for the numerical
flux to ensure positivity and to satisfy an entropy inequality.

Further details on this hydrostatic reconstruction and its motivation can be found in
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI:10.1137/15M1053074](https://dx.doi.org/10.1137/15M1053074)
"""
@inline function min_max_speed_chen_noelle(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquationsWetDry1D)
    # Get the velocity quantities
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Calculate the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)

    a_ll = sqrt(equations.gravity * h_ll)
    a_rr = sqrt(equations.gravity * h_rr)

    λ_min = min(v_ll - a_ll, v_rr - a_rr, zero(eltype(u_ll)))
    λ_max = max(v_ll + a_ll, v_rr + a_rr, zero(eltype(u_ll)))

    return λ_min, λ_max
end

@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.max_abs_speeds(u,
                                equations.basic_swe)
end

# Calculate estimates for minimum and maximum wave speeds for HLL-type fluxes
@inline function Trixi.min_max_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquationsWetDry1D)
    return Trixi.min_max_speed_naive(u_ll, u_rr, orientation,
                                     equations.basic_swe)
end

# More refined estimates for minimum and maximum wave speeds for HLL-type fluxes
@inline function Trixi.min_max_speed_davis(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquationsWetDry1D)
    return Trixi.min_max_speed_davis(u_ll, u_rr, orientation,
                                     equations.basic_swe)
end

@inline function Trixi.min_max_speed_einfeldt(u_ll, u_rr, orientation::Integer,
                                              equations::ShallowWaterEquationsWetDry1D)
    return Trixi.min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                        equations.basic_swe)
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.velocity(u,
                          equations.basic_swe)
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.cons2prim(u,
                           equations.basic_swe)
end

# Convert conservative variables to entropy
# Note, only the first two are the entropy variables, the third entry still
# just carries the bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.cons2entropy(u,
                              equations.basic_swe)
end

# Convert entropy variables to conservative
@inline function Trixi.entropy2cons(w, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.entropy2cons(w,
                              equations.basic_swe)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.prim2cons(prim,
                           equations.basic_swe)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterEquationsWetDry1D)
    return waterheight(u,
                       equations.basic_swe)
end

@inline function Trixi.pressure(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.pressure(u,
                          equations.basic_swe)
end

@inline function Trixi.waterheight_pressure(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.waterheight_pressure(u,
                                      equations.basic_swe)
end

# Entropy function for the shallow water equations is the total energy
@inline function Trixi.entropy(cons, equations::ShallowWaterEquationsWetDry1D)
    Trixi.energy_total(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline function Trixi.energy_total(cons, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.energy_total(cons,
                              equations.basic_swe)
end

# Calculate kinetic energy for a conservative state `cons`
@inline function Trixi.energy_kinetic(u, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.energy_kinetic(u,
                                equations.basic_swe)
end

# Calculate potential energy for a conservative state `cons`
@inline function Trixi.energy_internal(cons, equations::ShallowWaterEquationsWetDry1D)
    return Trixi.energy_total(cons, equations) - Trixi.energy_kinetic(cons, equations)
end

# Calculate the error for the "lake-at-rest" test case where H = h+b should
# be a constant value over time. Note, assumes there is a single reference
# water height `H0` with which to compare.
@inline function Trixi.lake_at_rest_error(u, equations::ShallowWaterEquationsWetDry1D)
    h, _, b = u

    # For well-balancedness testing with possible wet/dry regions the reference
    # water height `H0` accounts for the possibility that the bottom topography
    # can emerge out of the water as well as for the threshold offset to avoid
    # division by a "hard" zero water heights as well.
    H0_wet_dry = max(equations.H0, b + equations.threshold_limiter)

    return abs(H0_wet_dry - (h + b))
end
end # @muladd
