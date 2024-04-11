# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

####################################################################################################
# Include files with actual implementations for different systems of equations. 

include("shallow_water_wet_dry_1d.jl")
include("shallow_water_wet_dry_2d.jl")

include("shallow_water_two_layer_1d.jl")
include("shallow_water_two_layer_2d.jl")

abstract type AbstractShallowWaterMultiLayerEquations{NDIMS, NVARS, NLAYERS} <:
              Trixi.AbstractEquations{NDIMS, NVARS} end
include("shallow_water_multilayer_1d.jl")
include("shallow_water_multilayer_2d.jl")

"""
    eachlayer(equations::AbstractShallowWaterMultiLayerEquations)

Return an iterator over the indices that specify the location in relevant data structures
for the layers in `AbstractShallowWaterMultiLayerEquations`.
"""
@inline function eachlayer(equations::AbstractShallowWaterMultiLayerEquations)
    Base.OneTo(nlayers(equations))
end

"""
    nlayers(equations::AbstractShallowWaterMultiLayerEquations)

Retrieve the number of layers from an equation instance of the `AbstractShallowWaterMultiLayerEquations`.
"""
@inline function nlayers(::AbstractShallowWaterMultiLayerEquations{NDIMS, NVARS,
                                                                   NLAYERS}) where {
                                                                                    NDIMS,
                                                                                    NVARS,
                                                                                    NLAYERS
                                                                                    }
    NLAYERS
end
end # @muladd
