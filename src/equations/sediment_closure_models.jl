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

# Model dependent shear stress coefficient
@inline function shear_stress_coefficient(u, friction::ManningFriction)
    h, _, _, _ = u
    return friction.n^2 / h^(1 / 3)
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
end # @muladd
