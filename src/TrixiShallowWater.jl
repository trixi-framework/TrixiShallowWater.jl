module TrixiShallowWater

# While `using Trixi` makes all exported symbols available, in order to extend a method from the 
# `Trixi.jl` module, symbols need to be explicitly qualified with `Trixi.function_name`.
# For more information, see 
# https://github.com/trixi-framework/TrixiShallowWater.jl/pull/10#discussion_r1433720559
using Trixi
# Import additional symbols that are not exported by Trixi.jl
using Trixi: get_node_vars, set_node_vars!, waterheight
using MuladdMacro: @muladd
using StaticArrays: SVector, @SMatrix, MVector
using Static: True, False
using LinearAlgebra: norm

include("equations/equations.jl")
include("equations/numerical_fluxes.jl")
include("callbacks_stage/callbacks_stage.jl")
include("solvers/indicators.jl")

# Export types/functions that define the public API of TrixiShallowWater.jl
export ShallowWaterEquationsWetDry1D, ShallowWaterEquationsWetDry2D,
       ShallowWaterTwoLayerEquations1D, ShallowWaterTwoLayerEquations2D,
       ShallowWaterMultiLayerEquations1D, ShallowWaterMultiLayerEquations2D

export hydrostatic_reconstruction_chen_noelle, flux_nonconservative_chen_noelle,
       min_max_speed_chen_noelle, flux_hll_chen_noelle,
       flux_ersing_etal, flux_es_ersing_etal, hydrostatic_reconstruction_ersing_etal

export nlayers, eachlayer

export PositivityPreservingLimiterShallowWater

export IndicatorHennemannGassnerShallowWater

end
