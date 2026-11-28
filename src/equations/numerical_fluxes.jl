# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# This file contains general numerical fluxes that are not specific to certain equations

# An empty version of the `min_max_speed_chen_noelle` function is declared here
# in order to create a dimension agnostic version of `flux_hll_chen_noelle`.
# The full description of this wave speed estimate can be found in the docstrings
# for `min_max_speed_chen_noelle` in `shallow_water_wet_dry_1d.jl` or `shallow_water_wet_dry_2d.jl`.

function min_max_speed_chen_noelle end

"""
    flux_hll_chen_noelle = FluxHLL(min_max_speed_chen_noelle)

An instance of [`Trixi.FluxHLL`](@extref) specific to the shallow water equations that
uses the wave speed estimates from [`min_max_speed_chen_noelle`](@ref).
This HLL flux is guaranteed to have zero numerical mass flux out of a "dry" element,
maintain positivity of the water height, and satisfy an entropy inequality.

For complete details see Section 2.4 of the following reference
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI: 10.1137/15M1053074](https://doi.org/10.1137/15M1053074)
"""
const flux_hll_chen_noelle = FluxHLL(min_max_speed_chen_noelle)

"""
    FluxHydrostaticReconstruction(numerical_flux, hydrostatic_reconstruction)

Allow for some kind of hydrostatic reconstruction of the solution state prior to the
surface flux computation. This is a particular strategy to ensure that the method remains
well-balanced for the shallow water equations, see [`ShallowWaterEquations1D`](@ref)
or [`ShallowWaterEquations2D`](@ref).

For example, the hydrostatic reconstruction from Audusse et al. is implemented
in one and two spatial dimensions, see [`hydrostatic_reconstruction_audusse_etal`](@ref) or
the original paper
- Emmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004)
  A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
  [DOI: 10.1137/S1064827503431090](https://doi.org/10.1137/S1064827503431090)

Other hydrostatic reconstruction techniques are available, particularly to handle wet / dry
fronts. A good overview of the development and application of hydrostatic reconstruction can be found in
- Guoxian Chen and Sebastian Noelle
  A unified surface-gradient and hydrostatic reconstruction scheme for the shallow water equations (2021)
  [RWTH Aachen preprint](https://www.igpm.rwth-aachen.de/forschung/preprints/517)
- Andreas Buttinger-Kreuzhuber, Zsolt Horváth, Sebastian Noelle, Günter Blöschl and Jürgen Waser (2019)
  A fast second-order shallow water scheme on two-dimensional structured grids over abrupt topography
  [DOI: 10.1016/j.advwatres.2019.03.010](https://doi.org/10.1016/j.advwatres.2019.03.010)
"""
struct FluxHydrostaticReconstruction{NumericalFlux, HydrostaticReconstruction}
    numerical_flux::NumericalFlux
    hydrostatic_reconstruction::HydrostaticReconstruction
end

@inline function (numflux::FluxHydrostaticReconstruction)(u_ll, u_rr,
                                                          orientation_or_normal_direction,
                                                          equations::Trixi.AbstractEquations)
    @unpack numerical_flux, hydrostatic_reconstruction = numflux

    # Create the reconstructed left/right solution states in conservative form
    u_ll_star, u_rr_star = hydrostatic_reconstruction(u_ll, u_rr, equations)

    # Use the reconstructed states to compute the numerical surface flux
    return numerical_flux(u_ll_star, u_rr_star, orientation_or_normal_direction,
                          equations)
end

# Equations that store quantities that are constant in time, e.g., the bottom topography, in the
# state vector. Analogously to the specializations of `DissipationLocalLaxFriedrichs`, the
# global Lax-Friedrichs dissipation and flux splitting are specialized for these equations to
# avoid spurious dissipation in these quantities.
const EquationsWithStaticVariables = Union{ShallowWaterEquations1D,
                                           ShallowWaterEquations2D,
                                           ShallowWaterEquationsQuasi1D,
                                           ShallowWaterTwoLayerEquations1D,
                                           ShallowWaterTwoLayerEquations2D,
                                           ShallowWaterMultiLayerEquations1D,
                                           ShallowWaterMultiLayerEquations2D,
                                           ShallowWaterMomentEquations1D,
                                           ShallowWaterLinearizedMomentEquations1D,
                                           HyperbolicSainteMarieEquations1D}

# Number of trailing variables in the state vector that are constant in time
@inline nstaticvariables(::EquationsWithStaticVariables) = 1
# The bottom topography and the channel width
@inline nstaticvariables(::ShallowWaterEquationsQuasi1D) = 2

# Set the variables that are constant in time to zero
@inline function zero_static_variables(u, equations::EquationsWithStaticVariables)
    n_evolved = nvariables(equations) - nstaticvariables(equations)
    return SVector(ntuple(i -> i <= n_evolved ? u[i] : zero(eltype(u)),
                          Val(nvariables(equations))))
end

@inline function (dissipation::DissipationGlobalLaxFriedrichs)(u_ll, u_rr,
                                                               orientation::Integer,
                                                               equations::EquationsWithStaticVariables)
    @unpack λ = dissipation
    return -0.5f0 * λ * zero_static_variables(u_rr - u_ll, equations)
end

@inline function (dissipation::DissipationGlobalLaxFriedrichs)(u_ll, u_rr,
                                                               normal_direction::AbstractVector,
                                                               equations::EquationsWithStaticVariables)
    @unpack λ = dissipation
    return -0.5f0 * λ * norm(normal_direction) *
           zero_static_variables(u_rr - u_ll, equations)
end

# The flux vanishes in the variables that are constant in time, such that it suffices to set
# these variables to zero in the dissipative part of the splitting
@inline function (splitting::SplittingGlobalLaxFriedrichs)(u, ::Val{:plus},
                                                           orientation::Integer,
                                                           equations::EquationsWithStaticVariables)
    @unpack λ = splitting
    return 0.5f0 * (flux(u, orientation, equations) +
            λ * zero_static_variables(u, equations))
end

@inline function (splitting::SplittingGlobalLaxFriedrichs)(u, ::Val{:minus},
                                                           orientation::Integer,
                                                           equations::EquationsWithStaticVariables)
    @unpack λ = splitting
    return 0.5f0 * (flux(u, orientation, equations) -
            λ * zero_static_variables(u, equations))
end

@inline function (splitting::SplittingGlobalLaxFriedrichs)(u, ::Val{:plus},
                                                           normal_direction::AbstractVector,
                                                           equations::EquationsWithStaticVariables)
    @unpack λ = splitting
    return 0.5f0 * (flux(u, normal_direction, equations) +
            λ * norm(normal_direction) * zero_static_variables(u, equations))
end

@inline function (splitting::SplittingGlobalLaxFriedrichs)(u, ::Val{:minus},
                                                           normal_direction::AbstractVector,
                                                           equations::EquationsWithStaticVariables)
    @unpack λ = splitting
    return 0.5f0 * (flux(u, normal_direction, equations) -
            λ * norm(normal_direction) * zero_static_variables(u, equations))
end
end # @muladd
