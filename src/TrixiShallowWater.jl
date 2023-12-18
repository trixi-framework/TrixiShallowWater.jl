module TrixiShallowWater
# Decide between using Trixi: Trixi, import Trixi or using Trixi?
using Trixi: Trixi
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False

include("equations/equations.jl")

baz() = Trixi.examples_dir()

# export types/functions that define the public API of TrixiShallowWater.jl
export ShallowWaterEquationsWetDry1D
# TODO: These function are currently exported by Trixi.jl. Needs to be uncommented when removed from Trixi.jl
#export hydrostatic_reconstruction_chen_noelle, flux_nonconservative_chen_noelle, min_max_speed_chen_noelle

end
