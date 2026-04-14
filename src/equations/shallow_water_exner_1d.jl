# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Abstract type for the different bottom friction models
abstract type Friction{RealT} end

"""
    ManningFriction(; n)

Creates a Manning friction model for the bottom friction with Manning coefficient `n`.
The type is used to dispatch on the respective friction law through the `shear_stress_coefficient` 
when computing the `shear_stress`.
"""
struct ManningFriction{RealT} <: Friction{RealT}
    n::RealT
end

function ManningFriction(; n)
    ManningFriction(n)
end

# Abstract type for the different models to compute sediment discharge
abstract type SedimentModel{RealT} end

@doc raw"""
    ShieldsStressModel(; m_1, m_2, m_3, k_1, k_2, k_3, theta_c, d_s)

Create a Shields stress model to compute the sediment discharge `q_s` based on the generalized
formulation from equation (1.2) in the given reference.

The choice of the real constants `m_1`, `m_2`, `m_3`, `k_1`, `k_2`, and `k_3` creates
different models. For example, setting `m_1=0`, `m_2=1.5`, `m_3=0`, `k_1=8`, `k_2=1`, and `k_3=0`
yields the sedimentation model of Meyer-Peter and Müller as given in [`MeyerPeterMueller`](@ref) below.
The Shields stress represents the ratio of agitating and stabilizing forces in the sediment bed where
`theta_c` is the critical Shields stress for incipient motion and `d_s` is the mean diameter of
the sediment grain size.

- E.D. Fernández-Nieto, T.M. de Luna, G. Narbona-Reina and J. de Dieu Zabsonré (2017)\
  Formal deduction of the Saint-Venant–Exner model including arbitrarily sloping sediment beds and
  associated energy\
  [DOI: 10.1051/m2an/2016018](https://doi.org/10.1051/m2an/2016018)
"""
struct ShieldsStressModel{RealT} <: SedimentModel{RealT}
    m_1::RealT
    m_2::RealT
    m_3::RealT
    k_1::RealT
    k_2::RealT
    k_3::RealT
    theta_c::RealT    # critical shields stress
    d_s::RealT        # grain diameter
end

@doc raw"""
    GrassModel(; A_g, m_g=3)

Creates a Grass model to compute the sediment discharge `q_s` as
```math
q_s = A_g v^{m_g}
```
with the coefficients `A_g` and `m_g`. The constant `A_g` lies in the interval ``[0,1]``
and is a dimensional calibration constant that is usually measured experimentally.
It expresses the kind of interaction between the fluid and the sediment, the strength
of which increases as `A_g` approaches to 1. The factor `m_g` lies in the interval ``[1, 4]``.
Typically, one considers an odd integer value for `m_g` such that the sediment discharge
`q_s` can be differentiated and the model remains valid for all values of the velocity `v`.

An overview of different formulations to compute the sediment discharge can be found in:
- M.J. Castro Díaz, E.D. Fernández-Nieto, A.M. Ferreiro (2008)\
  Sediment transport models in Shallow Water equations and numerical approach by high order
  finite volume methods\
  [DOI:10.1016/j.compfluid.2007.07.017](https://doi.org/10.1016/j.compfluid.2007.07.017)
"""
struct GrassModel{RealT} <: SedimentModel{RealT}
    A_g::RealT
    m_g::RealT
end

function GrassModel(; A_g, m_g = 3)
    RealT = promote_type(typeof(A_g), typeof(m_g))
    return GrassModel(RealT(A_g), RealT(m_g))
end

@doc raw"""
    MeyerPeterMueller(; theta_c, d_s)

Creates a Meyer-Peter-Mueller model to compute the sediment discharge
`q_s` with the critical Shields stress `theta_c` and the grain diameter `d_s`.

An overview of different formulations to compute the sediment discharge can be found in:
- M.J. Castro Díaz, E.D. Fernández-Nieto, A.M. Ferreiro (2008)\
  Sediment transport models in Shallow Water equations and numerical approach by high order
  finite volume methods\
  [DOI:10.1016/j.compfluid.2007.07.017](https://doi.org/10.1016/j.compfluid.2007.07.017)
"""
function MeyerPeterMueller(; theta_c, d_s)
    RealT = promote_type(typeof(theta_c), typeof(d_s))
    return ShieldsStressModel(RealT(0.0), RealT(1.5), RealT(0.0), RealT(8.0),
                              RealT(1.0), RealT(0.0), RealT(theta_c), RealT(d_s))
end

@doc raw"""
    ShallowWaterExnerEquations1D(;gravity, H0 = 0.0,
                                 friction = ManningFriction(n = 0.0),
                                 sediment_model,
                                 porosity,
                                 rho_f, rho_s)
Formulation of the Shallow water-Exner equations in one space dimension that possesses a mathematical
entropy inequality. The equations are given by
```math
\begin{cases}
\partial_t h + \partial_x hv = 0, \\
\partial_t hv + \partial_x (hv^2) + gh\partial_x (h + h_b) + g\frac{1}{r}h_s\partial_x (rh + h_b) + \frac{\tau}{\rho_f} = 0,\\
\partial_t h_b + \partial_x q_s = 0,
\end{cases}
```
The unknown quantities are the water and sediment height ``h``, ``h_b`` and the velocity ``v``.
The sediment discharge ``q_s(h, hv)`` is determined by the `sediment_model` and is used to determine
the active sediment height ``h_s = q_s / v``.
Furthermore ``\tau`` denotes the shear stress at the water-sediment interface and is determined by
the `friction` model.
The gravitational acceleration is denoted by ``g``, and ``\rho_f`` and ``\rho_s`` are the fluid and sediment
densities, respectively. The density ratio is given by ``r = \rho_f / \rho_s``, where ``r`` lies between
``0 < r < 1`` as the fluid density ``\rho_f`` should be smaller than the sediment density ``\rho_s``.

The conservative variable water height ``h`` is measured from the sediment height ``h_b``, therefore
one also defines the total water height as ``H = h + h_b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water
height that is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

The entropy conservative formulation has been derived in the paper:
- E.D. Fernández-Nieto, T.M. de Luna, G. Narbona-Reina and J. de Dieu Zabsonré (2017)\
  Formal deduction of the Saint-Venant–Exner model including arbitrarily sloping sediment beds and
  associated energy\
  [DOI: 10.1051/m2an/2016018](https://doi.org/10.1051/m2an/2016018)
"""
struct ShallowWaterExnerEquations1D{RealT <: Real,
                                    FrictionT <: Friction{RealT},
                                    SedimentT <: SedimentModel{RealT}} <:
       Trixi.AbstractShallowWaterEquations{1, 3}
    # This structure ensures that friction and sediment models are concrete types
    # to prevent allocations
    gravity::RealT # gravitational acceleration
    H0::RealT      # constant "lake-at-rest" total water height
    friction::FrictionT
    sediment_model::SedimentT
    porosity_inv::RealT  # 1/(1-porosity)
    rho_f::RealT    # density of fluid
    rho_s::RealT    # density of sediment
    r::RealT       # density ratio
end

# Allow for flexibility to set the gravitational acceleration within an elixir depending on the
# application where `gravity=1.0` or `gravity=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases.
function ShallowWaterExnerEquations1D(; gravity, H0 = zero(gravity),
                                      friction = ManningFriction(n = 0.0),
                                      sediment_model,
                                      porosity, rho_f, rho_s)
    # Precompute common expressions for the porosity and density ratio
    porosity_inv = inv(1 - porosity)
    r = rho_f / rho_s
    return ShallowWaterExnerEquations1D(gravity, H0, friction, sediment_model,
                                        porosity_inv, rho_f, rho_s, r)
end

Trixi.have_nonconservative_terms(::ShallowWaterExnerEquations1D) = True()
Trixi.varnames(::typeof(cons2cons), ::ShallowWaterExnerEquations1D) = ("h", "hv", "h_b")
# Note, we use the total water height, H = h + h_b, as the first primitive variable for easier
# visualization and setting initial conditions
Trixi.varnames(::typeof(cons2prim), ::ShallowWaterExnerEquations1D) = ("H", "v", "h_b")

@doc raw"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, x, t, surface_flux_function,
                                  equations::ShallowWaterExnerEquations1D)

Create a boundary state by reflecting the normal velocity component and keep
the tangential velocity component unchanged. The boundary water height is taken from
the internal value.

For details see Section 9.2.5 of the book:
- Eleuterio F. Toro (2001)\
  Shock-Capturing Methods for Free-Surface Shallow Flows\
  1st edition\
  ISBN 0471987662
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation_or_normal,
                                                    direction,
                                                    x, t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterExnerEquations1D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

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

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterExnerEquations1D)

A smooth initial condition used for convergence tests in combination with
[`Trixi.source_terms_convergence_test`](@ref).
"""
@inline function Trixi.initial_condition_convergence_test(x, t,
                                                          equations::ShallowWaterExnerEquations1D)
    ω = sqrt(2) * pi

    h = 2 + cos(ω * x[1]) * cos(ω * t)
    v = 0.5f0
    h_b = 2 + sin(ω * x[1]) * cos(ω * t)

    return SVector(h, h * v, h_b)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterExnerEquations1D{T, S, GrassModel{T}}) where {T, S}

Source terms used for convergence tests in combination with [`Trixi.initial_condition_convergence_test`](@ref)
when using the the [`GrassModel`](@ref) model.

To use this source term the equations must be set to:
```julia
equations = ShallowWaterExnerEquations1D(gravity = 10.0, rho_f = 0.5,
                                            rho_s = 1.0, porosity = 0.5,
                                            friction = ManningFriction(n = 0.0),
                                            sediment_model = GrassModel(A_g = 0.01))
```
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterExnerEquations1D{T,
                                                                                             S,
                                                                                             GrassModel{T}}) where {
                                                                                                                    T,
                                                                                                                    S
                                                                                                                    }
    ω = sqrt(2) * pi
    A_g = equations.sediment_model.A_g

    h = -cos(x[1] * ω) * sin(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω
    hv = -0.5f0 * cos(x[1] * ω) * sin(t * ω) * ω -
         0.25f0 * sin(x[1] * ω) * cos(t * ω) * ω +
         10 * A_g *
         (cos(x[1] * ω) * cos(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω) +
         10 * (2 + cos(x[1] * ω) * cos(t * ω)) *
         (cos(x[1] * ω) * cos(t * ω) * ω - sin(x[1] * ω) * cos(t * ω) * ω)
    h_b = -sin(x[1] * ω) * sin(t * ω) * ω
    return SVector(h, hv, h_b)
end

""" 
    source_terms_convergence_test(u, x, t, equations::ShallowWaterExnerEquations1D{T, S, ShieldsStressModel{T}}) where {T, S}

Source terms used for convergence tests in combination with [`Trixi.initial_condition_convergence_test`](@ref)
when using the [`MeyerPeterMueller`](@ref) model.

To use this source term the equations must be set to:
```julia
equations = ShallowWaterExnerEquations1D(gravity = 10.0, rho_f = 0.5,
                                         rho_s = 1.0, porosity = 0.5,
                                         friction = ManningFriction(n = 0.01),
                                         sediment_model = MeyerPeterMueller(theta_c = 0.0,
                                                                            d_s = 1e-3))
```
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterExnerEquations1D{T,
                                                                                             S,
                                                                                             ShieldsStressModel{T}}) where {
                                                                                                                            T,
                                                                                                                            S
                                                                                                                            }
    ω = sqrt(2) * pi
    (; gravity, porosity_inv, rho_f, rho_s, r) = equations

    n = equations.friction.n

    # Constant expression from the MPM model
    c = sqrt(gravity * (1 / r - 1)) * 8 * porosity_inv *
        (rho_f / (rho_s - rho_f))^(3 / 2) * n^3

    h = -cos(x[1] * ω) * sin(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω

    hv = ((5 * c *
           (cos(x[1] * ω) * cos(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω)) /
          ((2 + cos(x[1] * ω) * cos(t * ω))^0.5) -
          0.5f0 * cos(x[1] * ω) * sin(t * ω) * ω -
          0.25f0 * sin(x[1] * ω) * cos(t * ω) * ω +
          10 * (2 + cos(x[1] * ω) * cos(t * ω)) *
          (cos(x[1] * ω) * cos(t * ω) * ω - sin(x[1] * ω) * cos(t * ω) * ω))

    h_b = ((0.5f0 * ((0.125f0 * c) / (2 + cos(x[1] * ω) * cos(t * ω))) * sin(x[1] * ω) *
            cos(t * ω) * ω) / ((2 + cos(x[1] * ω) * cos(t * ω))^0.5) -
           sin(x[1] * ω) * sin(t * ω) * ω)

    return SVector(h, hv, h_b)
end

"""
    source_term_bottom_friction(u, x, t, equations::ShallowWaterExnerEquations1D)

Source term that accounts for the bottom friction in the [ShallowWaterExnerEquations1D](@ref).
The actual friction law is determined through the friction model in `equations.friction`.
"""
@inline function source_term_bottom_friction(u, x, t,
                                             equations::ShallowWaterExnerEquations1D)
    return SVector(zero(eltype(u)),
                   -shear_stress(u, equations),
                   zero(eltype(u)))
end

# Calculate 1D flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::ShallowWaterExnerEquations1D)
    _, hv, _ = u

    v = velocity(u, equations)

    f1 = hv
    f2 = hv * v
    f3 = q_s(u, equations)

    return SVector(f1, f2, f3)
end

"""
    flux_nonconservative_ersing_etal(u_ll, u_rr, orientation, equations::ShallowWaterExnerEquations1D)

Non-symmetric path-conservative two-point flux discretizing the nonconservative terms of the
[`ShallowWaterExnerEquations1D`](@ref) which consists of the hydrostatic pressure of the fluid
layer and an additional pressure contribution from the sediment layer to obtain an entropy
inequality.

This non-conservative flux should be used together with [`flux_ersing_etal`](@ref) to create a
scheme that is entropy conservative and well-balanced.
"""
@inline function flux_nonconservative_ersing_etal(u_ll, u_rr,
                                                  orientation::Integer,
                                                  equations::ShallowWaterExnerEquations1D)
    # Pull the necessary left and right state information
    h_ll, _, h_b_ll = u_ll
    h_rr, _, h_b_rr = u_rr

    # Calculate jumps
    h_jump = h_rr - h_ll
    h_b_jump = h_b_rr - h_b_ll

    # Workaround to avoid division by zero, when computing the effective sediment height
    if abs(velocity(u_ll, equations)) < eps(typeof(h_ll))
        h_s_ll = 0
    else
        h_s_ll = q_s(u_ll, equations) / velocity(u_ll, equations)
    end

    z = zero(eltype(u_ll))

    f2 = equations.gravity * h_ll * (h_jump + h_b_jump)
    # Additional nonconservative term to obtain entropy conservative formulation
    f2 += (equations.gravity / equations.r * h_s_ll *
           (equations.r * h_jump + h_b_jump))

    return SVector(z, f2, z)
end

"""
    flux_ersing_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterMultiLayerEquations1D)

Entropy conservative split form, without the hydrostatic pressure. This flux should be used
together with the nonconservative [`flux_nonconservative_ersing_etal`](@ref) to create a scheme
that is entropy conservative and well-balanced.

To obtain an entropy stable formulation the `surface_flux` can be set as
`FluxPlusDissipation(flux_ersing_etal, DissipationLocalLaxFriedrichs()), flux_nonconservative_ersing_etal`.
"""
@inline function flux_ersing_etal(u_ll, u_rr, orientation::Integer,
                                  equations::ShallowWaterExnerEquations1D)
    # Unpack left and right state
    _, h_v_ll, _ = u_ll
    _, h_v_rr, _ = u_rr

    # Get the velocities on either side
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Average each factor of products in flux
    v_avg = 0.5f0 * (v_ll + v_rr)

    # Calculate fluxes depending on orientation
    f1 = 0.5f0 * (h_v_ll + h_v_rr)
    f2 = f1 * v_avg
    f3 = 0.5f0 * (q_s(u_ll, equations) + q_s(u_rr, equations))

    return SVector(f1, f2, f3)
end

"""
    dissipation_roe(u_ll, u_rr, orientation_or_normal_direction,
                                    equations::ShallowWaterExnerEquations1D)
Roe-type dissipation term for the [`ShallowWaterExnerEquations1D`](@ref) with an approximate Roe average
for the sediment discharge `q_s`.
"""
@inline function dissipation_roe(u_ll, u_rr, orientation_or_normal_direction,
                                 equations::ShallowWaterExnerEquations1D)
    r = equations.r
    g = equations.gravity
    z = zero(eltype(u_ll))

    # Get entropy variables and velocities
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Compute approximate Roe averages.
    # The actual Roe average for the sediment height `h_b` depends on the sediment and
    # friction model and an explicit formula is not always available.
    # Therefore we only use an approximation here.
    h_avg = 0.5f0 * (u_ll[1] + u_rr[1])
    v_avg = (sqrt(u_ll[1]) * v_ll + sqrt(u_rr[1]) * v_rr) /
            (sqrt(u_ll[1]) + sqrt(u_rr[1]))
    h_b_avg = 0.5f0 * (u_ll[3] + u_rr[3])

    # Workaround to avoid division by zero, when computing the effective sediment height
    if abs(v_avg) < eps(typeof(h_avg))
        h_s_avg = 0
    else
        h_s_avg = (q_s(SVector(h_avg, h_avg * v_avg, h_b_avg), equations) /
                   v_avg)
    end

    # Compute the eigenvalues using Cardano's formula
    λ1, λ2, λ3 = eigvals_cardano(SVector(h_avg, h_avg * v_avg, h_b_avg), equations)

    # Precompute some common expressions
    c1 = g * (h_avg + h_s_avg)
    c2 = g * (h_avg + h_s_avg / r)

    # Eigenvector matrix
    r31 = ((v_avg - λ1)^2 - c1) / c2
    r32 = ((v_avg - λ2)^2 - c1) / c2
    r33 = ((v_avg - λ3)^2 - c1) / c2
    R = @SMatrix [[1 1 1]; [λ1 λ2 λ3]; [r31 r32 r33]]

    # Inverse eigenvector matrix
    d1 = (λ1 - λ2) * (λ1 - λ3)
    d2 = (λ2 - λ1) * (λ2 - λ3)
    d3 = (λ3 - λ2) * (λ3 - λ1)
    R_inv = @SMatrix [(c1 - v_avg^2 + λ2 * λ3)/d1 (2 * v_avg - λ2 - λ3)/d1 c2/d1;
                      (c1 - v_avg^2 + λ1 * λ3)/d2 (2 * v_avg - λ1 - λ3)/d2 c2/d2;
                      (c1 - v_avg^2 + λ1 * λ2)/d3 (2 * v_avg - λ2 - λ1)/d3 c2/d3]

    # Eigenvalue absolute value matrix
    Λ_abs = @SMatrix [abs(λ1) z z; z abs(λ2) z; z z abs(λ3)]

    # Compute dissipation
    diss = SVector(-0.5f0 * R * Λ_abs * R_inv * (u_rr - u_ll))

    return SVector(diss[1], diss[2], diss[3])
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterExnerEquations1D)
    return max(maximum(abs, eigvals_cardano(u_rr, equations)),
               maximum(abs, eigvals_cardano(u_ll, equations)))
end

@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterExnerEquations1D)
    return maximum(abs, eigvals_cardano(u, equations))
end

#Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterExnerEquations1D)
    h, hv, _ = u

    v = hv / h

    return v
end

# Compute the sediment discharge for Shields stress models
@inline function q_s(u,
                     equations::ShallowWaterExnerEquations1D{T, S,
                                                             ShieldsStressModel{T}}) where {
                                                                                            T,
                                                                                            S
                                                                                            }
    (; gravity, rho_f, rho_s, porosity_inv) = equations
    (; m_1, m_2, m_3, k_1, k_2, k_3, theta_c, d_s) = equations.sediment_model

    theta = rho_f * abs(shear_stress(u, equations)) / (gravity * (rho_s - rho_f) * d_s)  # Shields stress

    Q = d_s * sqrt(gravity * (rho_s / rho_f - 1) * d_s) # Characteristic discharge

    return (porosity_inv * Q * sign(theta) * k_1 * theta^m_1 *
            (max(theta - k_2 * theta_c, 0))^m_2 *
            (max(sqrt(theta) - k_3 * sqrt(theta_c), 0))^m_3)
end

# Compute the sediment discharge for the Grass model
@inline function q_s(u,
                     equations::ShallowWaterExnerEquations1D{T, S, GrassModel{T}}) where {
                                                                                          T,
                                                                                          S
                                                                                          }
    (; porosity_inv, sediment_model) = equations
    return porosity_inv * sediment_model.A_g * velocity(u, equations)^sediment_model.m_g
end

# Shear stress formulation using a coefficient to take into account different friction models
@inline function shear_stress(u, equations::ShallowWaterExnerEquations1D)
    v = velocity(u, equations)
    return equations.gravity * shear_stress_coefficient(u, equations.friction) * v *
           abs(v)
end

# Model dependent shear stress coefficient
@inline function shear_stress_coefficient(u, friction::ManningFriction)
    h, _, _ = u
    return friction.n^2 / h^(1 / 3)
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterExnerEquations1D)
    h, _, h_b = u

    H = h + h_b
    v = velocity(u, equations)

    return SVector(H, v, h_b)
end

# Convert conservative variables to entropy variables
@inline function Trixi.cons2entropy(u, equations::ShallowWaterExnerEquations1D)
    h, _, h_b = u
    (; gravity, r) = equations

    v = velocity(u, equations)

    w1 = r * (gravity * (h + h_b) - 0.5f0 * v^2)
    w2 = r * v
    w3 = gravity * (r * h + h_b)

    return SVector(w1, w2, w3)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterExnerEquations1D)
    H, v, h_b = prim

    h = H - h_b
    hv = h * v

    return SVector(h, hv, h_b)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterExnerEquations1D)
    return u[1]
end

# Indicator variable used for the shock capturing in `IndicatorHennemannGassnerShallowWater`
@inline function water_sediment_height(u, equations::ShallowWaterExnerEquations1D)
    return equations.gravity * u[1] * u[3]
end

# Mathematical entropy function
@inline function Trixi.entropy(u, equations::ShallowWaterExnerEquations1D)
    h, _, h_b = u
    v = velocity(u, equations)
    (; gravity, r) = equations

    return 0.5f0 * r * h * v^2 + 0.5f0 * gravity * (r * h^2 + h_b^2) +
           r * gravity * h * h_b
end

# Calculate the error for the "lake-at-rest" test case where H = h + h_b should
# be a constant value over time.
@inline function Trixi.lake_at_rest_error(u, equations::ShallowWaterExnerEquations1D)
    h, _, h_b = u
    return abs(equations.H0 - (h + h_b))
end

# Trigonometric version of Cardano's method for the roots of a cubic polynomial
# in order to compute the eigenvalues of the [`ShallowWaterExnerEquations1D`[(@ref)].
# Note, assumes only real roots.
@inline function eigvals_cardano(u, equations::ShallowWaterExnerEquations1D)
    h = waterheight(u, equations)
    v = velocity(u, equations)
    g = equations.gravity
    r = equations.r

    # Workaround to avoid division by zero, when computing the effective sediment height
    if abs(v) > eps(typeof(h))
        h_s = q_s(u, equations) / v
        # Compute gradients of q_s using automatic differentiation.
        # Introduces a closure to make q_s a function of u only. This is necessary since the
        # gradient function only accepts functions of one variable.
        dq_s_dh, dq_s_dhv, _ = Trixi.ForwardDiff.gradient(u -> q_s(u, equations), u)
    else
        h_s = 0
        dq_s_dh = 0
        dq_s_dhv = 0
    end

    # Coefficient for the original cubic equation ax^3 + bx^2 + cx + d
    a = -1
    b = 2 * v
    c = g * (h + 1 / r * h_s) * dq_s_dhv + g * (h + h_s) - v^2
    d = g * (h + 1 / r * h_s) * dq_s_dh

    # Coefficient of the depressed cubic equation t^3 + pt + q = 0
    p = (3 * a * c - b^2) / (3 * a^2)
    q = (2 * b^3 - 9 * a * b * c + 27 * a^2 * d) / (27 * a^3)

    # Roots of the original cubic equation
    λ1 = -b / (3 * a) +
         2 * sqrt(-p / 3) * cos(1 / 3 * acos(3 * q / (2 * p) * sqrt(-3 / p)))
    λ2 = -b / (3 * a) +
         2 * sqrt(-p / 3) *
         cos(1 / 3 * acos(3 * q / (2 * p) * sqrt(-3 / p)) - 2 * π * 1 / 3)
    λ3 = -b / (3 * a) +
         2 * sqrt(-p / 3) *
         cos(1 / 3 * acos(3 * q / (2 * p) * sqrt(-3 / p)) - 2 * π * 2 / 3)

    return SVector(λ1, λ2, λ3)
end
end # @muladd
