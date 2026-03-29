@muladd begin
#! format: noindent

@doc raw"""
    HyperbolicSainteMarieEquations1D(; bathymetry_type = bathymetry_mild_slope,
                                       gravity,
                                       eta0 = 0.0,
                                       h0,
                                       alpha = 3.0)

Hyperbolic approximation of the Sainte-Marie system
[`SainteMarieEquations1D`](@ref) in one spatial dimension
(with parameter ``\gamma = 2`` compared to the original literature)
derived by Escalante, Dumbser and Castro (2019).
The equations are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + (h p)_x &= -(g h + 2 p) b_x,\\
  h w_t + h v w_x &= 2 p,\\
  h p_t + h v p_x + h c^2 \bigl(v_x + (w - v b_x) / (h / 2)\bigr) &= 0.
\end{aligned}
```

The additional quantity ``H_0`` is also available to store a reference value for the total water height that
is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

The unknown quantities of the Sainte-Marie equations 
are the water height ``h`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the (possibly) variable bottom topography
(bathymetry) ``b(x)``. Conservative variable water height ``h`` is measured from the bottom topography ``b``, therefore one
also defines the total water height as ``H = h + b``.
There are two auxiliary variables:``w \approx -h v_x / 2 + v b_x`` and the non-hydrostatic pressure ``p``.
In the formal limit ``c^2 \to \infty``, the hyperbolic approximation recovers the original Sainte-Marie system.
Escalante, Dumbser and Castro (2019) choose the hyperbolization parameter as ``c = \alpha \sqrt{g h_0}`` for some background water height ``h_0``.
Thus, the hyperbolization parameter ``c^2`` is set by the keyword arguments `alpha` (``\alpha``), `gravity` (``g``), and `h0` (``h_0``).
The larger the value of ``\alpha``, the better the approximation of the original system, but also the stiffer the system.
Typically, ``\alpha`` should be chosen larger than 1 to get a good approximation of the original system, but the system also becomes stiffer for larger values of ``\alpha``.
Escalante, Dumbser and Castro (2019) often use the value ``\alpha = 3`` in their numerical experiments, which is also the default value in `HyperbolicSainteMarieEquations1D`.
For the special case ``\alpha = 0`` and with initial conditions ``p = 0 = w``, the system reduces to the classical shallow water equations.

References for the Sainte-Marie system and its hyperbolization can be found in
- Sainte-Marie (2011)
  Vertically averaged models for the free surface non-hydrostatic
  Euler system: derivation and kinetic interpretation
  [DOI: 10.1142/S0218202511005118](https://doi.org/10.1142/S0218202511005118)
- Bristeau, Mangeney, Sainte-Marie and Seguin (2015)
  An energy-consistent depth-averaged Euler system:
  derivation and properties
  [DOI: 10.3934/dcdsb.2015.20.961](https://doi.org/10.3934/dcdsb.2015.20.961)
- Aïssiouene, Bristeau, Godlewski, Mangeney, Parés Madroñal and Sainte-Marie (2020)
  A two-dimensional method for a family of dispersive shallow water models
  [DOI: 10.5802/smai-jcm.66](https://doi.org/10.5802/smai-jcm.66)
- Escalante, Dumbser and Castro (2019)
  An efficient hyperbolic relaxation system for dispersive non-hydrostatic
  water waves and its solution with high order discontinuous Galerkin schemes
  [DOI: 10.1016/j.jcp.2019.05.035](https://doi.org/10.1016/j.jcp.2019.05.035)
"""

struct HyperbolicSainteMarieEquations1D{RealT <: Real} <:
       Trixi.AbstractShallowWaterEquations{1, 5}
    gravity::RealT
    H0::RealT
    celerity::RealT
    threshold_limiter::RealT
end

function HyperbolicSainteMarieEquations1D(; gravity, H0 = zero(gravity),
                                          b0 = one(gravity), alpha = 3,
                                          threshold_limiter = nothing)
    T = promote_type(typeof(gravity), typeof(H0), typeof(b0), typeof(alpha))
    if threshold_limiter === nothing
        threshold_limiter = 500 * eps(T)
    end
    celerity = alpha * sqrt(gravity * b0)

    HyperbolicSainteMarieEquations1D(gravity, H0, celerity, threshold_limiter)
end

Trixi.have_nonconservative_terms(::HyperbolicSainteMarieEquations1D) = Trixi.True()

function Trixi.varnames(::typeof(cons2cons), ::HyperbolicSainteMarieEquations1D)
    ("h", "h_v", "h_w", "h_p", "b")
end
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
function Trixi.varnames(::typeof(cons2prim), ::HyperbolicSainteMarieEquations1D)
    ("H", "v", "w", "p", "b")
end

"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, x, t, surface_flux_function,
                                  equations::HyperbolicSainteMarieEquations1D)

Create a boundary state by reflecting the normal momentum component and keep
the tangential momentum component unchanged.

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
                                                    equations::HyperbolicSainteMarieEquations1D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # This can not be dispatched, when Flux Hydrostactic reconstruction is used
    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         -u_inner[2],
                         u_inner[3], u_inner[4], u_inner[5])

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

@inline function source_term_hyperbolic_sainte_marie(u, x, t,
                                                     equations::HyperbolicSainteMarieEquations1D)
    h, h_v, h_w, h_p, b = u

    p = h_p / h
    w = h_w / h

    du1 = zero(eltype(u))
    du2 = zero(eltype(u))
    du3 = 2 * p
    du4 = -2 * equations.celerity^2 * w
    du5 = zero(eltype(u))

    return SVector(du1, du2, du3, du4, du5)
end

"""
	flux_conservative_artiano_ranocha(u_ll, u_rr, normal_direction::AbstractVector, equations::HyperbolicSainteMarieEquations1D)

Total energy conserving and well-balanced two-point flux by
-  Marco Artiano, Hendrik Ranocha (2026)
   On Affordable High-Order Entropy-Conservative/Stable and 
   Well-Balanced Methods for Nonconservative Hyperbolic Systems
   [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
struct flux_conservative_artiano_ranocha{RealT <: Real}
    alpha_1::RealT
    alpha_2::RealT
    alpha_3::RealT
end

@inline function (flux_ec::flux_conservative_artiano_ranocha)(u_ll, u_rr,
                                                              orientation::Integer,
                                                              equations::HyperbolicSainteMarieEquations1D)
    alpha_1 = flux_ec.alpha_1
    alpha_2 = flux_ec.alpha_2
    alpha_3 = flux_ec.alpha_3

    # Pull the necessary left and right state information
    h_ll, h_v_ll, h_w_ll, h_p_ll, b_ll = u_ll
    h_rr, h_v_rr, h_w_rr, h_p_rr, b_rr = u_rr

    v_ll = h_v_ll / h_ll
    w_ll = h_w_ll / h_ll
    p_ll = h_p_ll / h_ll

    v_rr = h_v_rr / h_rr
    w_rr = h_w_rr / h_rr
    p_rr = h_p_rr / h_rr

    v_avg = 0.5f0 * (v_ll + v_rr)
    w_avg = 0.5f0 * (w_ll + w_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    h_avg = 0.5f0 * (h_ll + h_rr)
    h_v_avg = 0.5f0 * (h_v_ll + h_v_rr)
    h_p_avg = 0.5f0 * (h_p_ll + h_p_rr)
    h2_avg = 0.5f0 * (h_ll^2 + h_rr^2)

    f1 = alpha_1 * h_avg * v_avg + (1 - alpha_1) * h_v_avg
    pressure_terms = equations.gravity * (1 - alpha_1) * h_avg^2 +
                     equations.gravity * (alpha_1 - 0.5f0) * h2_avg +
                     alpha_3 * h_avg * p_avg + (1 - alpha_3) * h_p_avg
    f2 = f1 * v_avg + pressure_terms
    f3 = f1 * w_avg
    f4 = f1 * p_avg

    return SVector(f1,
                   f2,
                   f3, f4, zero(eltype(u_ll)))
end

"""
	flux_nonconservative_artiano_ranocha(u_ll, u_rr, normal_direction::AbstractVector, equations::HyperbolicSainteMarieEquations1D)

Total energy conserving and well-balanced two-point flux by
-  Marco Artiano, Hendrik Ranocha (2026)
   On Affordable High-Order Entropy-Conservative/Stable and 
   Well-Balanced Methods for Nonconservative Hyperbolic Systems
   [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
struct flux_nonconservative_artiano_ranocha{RealT <: Real}
    alpha_1::RealT
    alpha_2::RealT
    alpha_3::RealT
end

@inline function (flux_ec::flux_nonconservative_artiano_ranocha)(u_ll, u_rr,
                                                                 orientation::Integer,
                                                                 equations::HyperbolicSainteMarieEquations1D)
    alpha_1 = flux_ec.alpha_1
    alpha_2 = flux_ec.alpha_2
    alpha_3 = flux_ec.alpha_3
    # Pull the necessary left and right state information
    h_ll, h_v_ll, h_w_ll, h_p_ll, b_ll = u_ll
    h_rr, h_v_rr, h_w_rr, h_p_rr, b_rr = u_rr

    v_ll = h_v_ll / h_ll
    p_ll = h_p_ll / h_ll
    b_jump = b_rr - b_ll
    v_rr = h_v_rr / h_rr
    p_rr = h_p_rr / h_rr

    v_avg = 0.5f0 * (v_ll + v_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    h_avg = 0.5f0 * (h_ll + h_rr)

    v_jump = v_rr - v_ll

    f1 = zero(eltype(u_ll))
    f2_fluxdiff = alpha_1 * equations.gravity * h_avg * b_jump +
                  alpha_2 * 2 * p_avg * b_jump
    f2_pointwise = (1 - alpha_1) * equations.gravity * h_ll * b_jump +
                   (1 - alpha_2) * 2 * p_ll * b_jump
    f2 = f2_fluxdiff + f2_pointwise
    f3 = zero(eltype(u_ll))

    f4_fluxdiff = alpha_3 * equations.celerity^2 * h_avg * v_jump -
                  2 * alpha_2 * equations.celerity^2 * v_avg * b_jump
    f4_pointwise = (1 - alpha_3) * equations.celerity^2 * h_ll * v_jump -
                   2 * (1 - alpha_2) * equations.celerity^2 * v_ll * b_jump
    f4 = f4_fluxdiff + f4_pointwise

    return SVector(f1,
                   f2,
                   f3,
                   f4,
                   zero(eltype(u_ll)))
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function Trixi.max_abs_speed(u_ll, u_rr, orientation::Integer,
                                     equations::HyperbolicSainteMarieEquations1D)
    # Get the velocity quantities
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Calculate the wave celerity on the left and right
    h_ll = u_ll[1]
    h_rr = u_rr[1]
    p_ll = u_ll[4] / h_ll
    p_rr = u_rr[4] / h_rr
    c_ll = sqrt(equations.gravity * h_ll + p_ll + equations.celerity^2)
    c_rr = sqrt(equations.gravity * h_rr + p_rr + equations.celerity^2)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

@inline function Trixi.max_abs_speeds(u, equations::HyperbolicSainteMarieEquations1D)
    h = u[1]
    v = velocity(u, equations)
    p = u[4] / h
    c = sqrt(equations.gravity * h + p + equations.celerity^2)
    return (abs(v) + c,)
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::HyperbolicSainteMarieEquations1D)
    h, h_v, _ = u

    v = h_v / h

    return v
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::HyperbolicSainteMarieEquations1D)
    h, h_v, h_w, h_p, b = u

    H = h + b
    v = h_v / h
    w = h_w / h
    p = h_p / h
    return SVector(H, v, w, p, b)
end

# Convert conservative variables to entropy
# Note, only the first two are the entropy variables, the third entry still
# just carries the bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::HyperbolicSainteMarieEquations1D)
    h, h_v, h_w, h_p, b = u

    v = h_v / h
    w = h_w / h
    p = h_p / h

    w1 = equations.gravity * (h + b) - 0.5f0 * (v^2 + w^2 + p^2 / equations.celerity^2)
    w2 = v
    w3 = w
    w4 = p / equations.celerity^2
    return SVector(w1, w2, w3, w4, zero(eltype(u)))
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::HyperbolicSainteMarieEquations1D)
    H, v, w, p, b = prim

    h = H - b
    h_v = h * v
    h_w = h * w
    h_p = h * p

    return SVector(h, h_v, h_w, h_p, b)
end

# Entropy function for the shallow water equations is the total energy
@inline function Trixi.entropy(cons, equations::HyperbolicSainteMarieEquations1D)
    energy_total(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline function Trixi.energy_total(cons, equations::HyperbolicSainteMarieEquations1D)
    h, h_v, h_w, h_p, b = cons

    e = (h_v^2 + h_w^2) / (2 * h) + h_p^2 / (2 * equations.celerity^2 * h) +
        0.5f0 * equations.gravity * h^2 + equations.gravity * h * b
    return e
end

@inline function Trixi.lake_at_rest_error(u,
                                          equations::HyperbolicSainteMarieEquations1D)
    h, _, _, _, b = u
    H0_wet_dry = max(equations.H0, b + equations.threshold_limiter)

    return abs(H0_wet_dry - (h + b))
end
end # @muladd
