module TrixiShallowWater
# We decided to import only Trixi.jl and qualify symbols explicitly with e.g. `Trixi.function_name`.
# For more information, see 
# https://github.com/trixi-framework/TrixiShallowWater.jl/pull/10#discussion_r1433720559
using Trixi
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False

include("equations/equations.jl")

# export types/functions that define the public API of TrixiShallowWater.jl
export ShallowWaterEquationsWetDry1D, ShallowWaterEquationsWetDry2D

export hydrostatic_reconstruction_chen_noelle, flux_nonconservative_chen_noelle, min_max_speed_chen_noelle

end
