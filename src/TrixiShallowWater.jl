module TrixiShallowWater
# TODO: rewrite this
# While we do using Trixi.jl to extend a method from Trixi.jl symbols need to be qualified explicitly 
# e.g. `Trixi.function_name`.
# For more information, see 
# https://github.com/trixi-framework/TrixiShallowWater.jl/pull/10#discussion_r1433720559
using Trixi
# Import additional symbols that are not exported by Trixi
import Trixi: get_node_vars, set_node_vars!
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False
using LinearAlgebra: norm

include("equations/equations.jl")
include("equations/numerical_fluxes.jl")
include("callbacks_stage/callbacks_stage.jl")
include("solvers/indicators.jl")

# export types/functions that define the public API of TrixiShallowWater.jl
export ShallowWaterEquationsWetDry1D, ShallowWaterEquationsWetDry2D

export hydrostatic_reconstruction_chen_noelle, flux_nonconservative_chen_noelle,
       min_max_speed_chen_noelle,
       flux_hll_chen_noelle

export PositivityPreservingLimiterShallowWater
export IndicatorHennemannGassnerShallowWater

end
