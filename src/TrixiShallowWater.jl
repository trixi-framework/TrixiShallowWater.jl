module TrixiShallowWater

# While `using Trixi` makes all exported symbols available, in order to extend a method from the
# `Trixi.jl` module, symbols need to be explicitly qualified with `Trixi.function_name`.
# For more information, see
# https://github.com/trixi-framework/TrixiShallowWater.jl/pull/10#discussion_r1433720559
using Trixi
# Import additional symbols that are not exported by Trixi.jl
using Trixi: get_node_vars, set_node_vars!
using MuladdMacro: @muladd
using StaticArrays: SVector, @SMatrix, MVector
using Static: True, False
using LinearAlgebra: norm
using Roots: Order2, solve, ZeroProblem

include("equations/equations.jl")
include("equations/numerical_fluxes.jl")
include("solvers/indicators.jl")
include("solvers/dgsem_p4est/containers.jl")
include("solvers/dgsem_p4est/dg_2d.jl")
include("callbacks_stage/callbacks_stage.jl")
include("callbacks_step/callbacks_step.jl")

# Export types/functions that define the public API of TrixiShallowWater.jl
export ShallowWaterEquations1D, ShallowWaterEquations2D,
       ShallowWaterExnerEquations1D,
       ShallowWaterTwoLayerEquations1D, ShallowWaterTwoLayerEquations2D,
       ShallowWaterMultiLayerEquations1D, ShallowWaterMultiLayerEquations2D,
       ShallowWaterEquationsQuasi1D

export hydrostatic_reconstruction_chen_noelle, flux_nonconservative_chen_noelle,
       min_max_speed_chen_noelle, flux_hll_chen_noelle,
       flux_ersing_etal, flux_nonconservative_ersing_etal,
       flux_es_ersing_etal, hydrostatic_reconstruction_ersing_etal,
       flux_nonconservative_audusse_etal, hydrostatic_reconstruction_audusse_etal,
       FluxHydrostaticReconstruction

export ManningFriction, MeyerPeterMueller, GrassModel, ShieldsStressModel,
       dissipation_roe, water_sediment_height, source_term_bottom_friction

export BoundaryConditionWaterHeight, BoundaryConditionMomentum

export nlayers, eachlayer

export PositivityPreservingLimiterShallowWater

export IndicatorHennemannGassnerShallowWater

end
