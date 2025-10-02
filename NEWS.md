# Changelog

TrixiShallowWater.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes in the v0.2 lifecycle

#### Added
- Introduce node-wise limiting functionality for the `ShallowWaterMultiLayerEquations2D`. ([#111])

#### Changed
- Velocity desingularization procedure has been moved into a distinct `VelocityDesingularization` 
  stage callback. ([#111])

#### Deprecated

#### Removed

## Changes when updating to v0.2 from v0.1.x

#### Added
- New equation `ShallowWaterEquationsQuasi1D` and functions `FluxHydrostaticReconstruction`, 
  `flux_nonconservative_audusse_etal`, and `hydrostatic_reconstruction_audusse_etal` are now available
  through TrixiShallowWater.jl instead of Trixi.jl. ([#96])
  
#### Changed
- `ShallowWaterEquationsWetDry` have been renamed to `ShallowWaterEquations`. The source code
  for these equations is now implemented directly in TrixiShallowWater.jl ([#96]).

#### Deprecated

#### Removed

## Changes in the v0.1 lifecycle

#### Added

- Experimental support for well-balanced mortars together with AMR [#45]
- New boundary conditions `BoundaryConditionWaterHeight` and `BoundaryConditionMomentum` now 
  available for the `ShallowWaterEquationsWetDry` to impose either the water height or the momentum 
  at the boundary. [#91]

#### Changed

#### Deprecated

#### Removed
