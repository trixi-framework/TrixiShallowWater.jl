# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# TODO: Right now this is only implemented for the SWE. Should we move it to the equation specific folder?
# Should we have BC(h_boundary, equations) instead of BC(h_boundary)?
"""
    BoundaryConditionWaterHeight(h_boundary)

Create a boundary condition that uses `h_boundary` to specify a fixed water height at the 
boundary and extrapolates the velocity from the incoming Riemann invariant.

The external water height `h_boundary` can be specified as a constant value or as a function of time, e.g.
```julia
   BoundaryConditionWaterHeight(2.0)
   BoundaryConditionWaterHeight(t -> 2 + cos(t))
```

More details can be found in the paper:
- Lixiang Song, Jianzhong Zhou, Jun Guo, Qiang Zou, Yi Liu (2011)
  A robust well-balanced finite volume model for shallow water flows
  with wetting and drying over irregular terrain
  [doi: 10.1016/j.advwatres.2011.04.017](https://doi.org/10.1016/j.advwatres.2011.04.017)

!!! warning "Experimental code"
    This is an experimental feature and can change any time.
"""
struct BoundaryConditionWaterHeight{F}
    h_boundary::F
end

# If `h_boundary` is provided as a constant create a function `h_boundary(t) = h_boundary`
function BoundaryConditionWaterHeight(h_boundary::Real)
    return BoundaryConditionWaterHeight(t -> h_boundary)
end

"""
    BoundaryConditionMomentum(hv1_boundary, hv2_boundary)

Create a boundary condition that sets a fixed momentum in x- and y- direction, `hv1_boundary`
and `hv2_boundary`, at the boundary and extrapolates the water height `h_boundary` from the incoming
Riemann invariant.

The external momentum can be specified as a constant value or as a function of time, e.g.
```julia
   BoundaryConditionMomentum(2.0, 1.0)
   BoundaryConditionMomentum(t -> 2 + cos(t), t -> 1 + sin(t))
```

More details can be found in the paper:
- Lixiang Song, Jianzhong Zhou, Jun Guo, Qiang Zou, Yi Liu (2011)
  A robust well-balanced finite volume model for shallow water flows
  with wetting and drying over irregular terrain
  [doi: 10.1016/j.advwatres.2011.04.017](https://doi.org/10.1016/j.advwatres.2011.04.017)

!!! warning "Experimental code"
    This is an experimental feature and can change any time.
"""
struct BoundaryConditionMomentum{F1, F2}
    hv1_boundary::F1
    hv2_boundary::F2
end

# If `hv_boundary` is provided as a constant create a function `hv_boundary(t) = hv_boundary``
function BoundaryConditionMomentum(hv1_boundary::Real, hv2_boundary::Real)
    return BoundaryConditionMomentum((t -> hv1_boundary), (t -> hv2_boundary))
end

####################################################################################################
# Include files with actual implementations for different systems of equations. 

include("shallow_water_wet_dry_1d.jl")
include("shallow_water_wet_dry_2d.jl")

include("shallow_water_exner_1d.jl")

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

# TODO: Add suitable default thresholds for Float32
# Provide default thresholds dependent on number format (Currently default thresholds are only provided
# for Float64)
default_threshold_partially_wet(::Type{Float64}) = 1e-4
default_threshold_partially_wet(catchall) = throw(ArgumentError("threshold_partially_wet must be provided for non-Float64 types"))

default_threshold_desingularization(::Type{Float64}) = 1e-10
default_threshold_desingularization(catchall) = throw(ArgumentError("threshold_desingularization must be provided for non-Float64 types"))
end # @muladd
