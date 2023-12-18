module TrixiShallowWater
# Decide between using Trixi: Trixi, import Trixi or using Trixi?
using Trixi: Trixi
using MuladdMacro: @muladd
using StaticArrays: SVector
using Static: True, False

include("equations/equations.jl")

# export types/functions that define the public API of TrixiShallowWater.jl
export ShallowWaterEquationsWetDry1D
export hydrostatic_reconstruction_chen_noelle, flux_nonconservative_chen_noelle

end
