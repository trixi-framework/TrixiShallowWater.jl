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

# Additional version of `FluxHydrostaticReconstruction` to add support for nonconservative fluxes on
# unstructured meshes. These can depend on both the contravariant vectors (normal direction) at 
# the current node and the averaged ones.
@inline function (numflux::Trixi.FluxHydrostaticReconstruction)(u_ll, u_rr,
                                                                normal_direction_ll,
                                                                normal_direction_average,
                                                                equations::Trixi.AbstractEquations)
    @unpack numerical_flux, hydrostatic_reconstruction = numflux

    # Create the reconstructed left/right solution states in conservative form
    u_ll_star, u_rr_star = hydrostatic_reconstruction(u_ll, u_rr, equations)

    # Use the reconstructed states to compute the nonconservative surface flux
    return numerical_flux(u_ll_star, u_rr_star, normal_direction_ll,
                          normal_direction_average,
                          equations)
end
end
