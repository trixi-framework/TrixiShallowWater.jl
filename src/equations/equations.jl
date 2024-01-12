# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

####################################################################################################
# Include files with actual implementations for different systems of equations. 

include("shallow_water_wet_dry_1d.jl")
end # @muladd
