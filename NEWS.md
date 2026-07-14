# Changelog

TrixiShallowWater.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes in the v0.2 lifecycle

#### Added
- Introduce node-wise limiting functionality for the `ShallowWaterMultiLayerEquations2D`. ([#111])
- New equation types `ShallowWaterMomentEquations1D` and `ShallowWaterLinearizedMomentEquations1D` have been added. ([#128])
- `ShallowWaterExner` extended to 2D on `TreeMesh` ([#150]) and curvilinear meshes ([#159]).
- Experimental support for rainfall & soil infiltration source terms for `ShallowWaterEquations1D` and `ShallowWaterEquations2D`. ([#158])
- New variants of `limiter_shallow_water!` needed to ensure positivity preservation after coarsening and refinement steps in AMR via the keyword argument `limiter!` in `AMRCallback` ([#164]). For details on this callback see ([#2396](https://github.com/trixi-framework/Trixi.jl/pull/2396)) in Trixi.jl. This made the specialized `refine!` and `coarsen!` introduced in [#97] obsolete, which is why they were removed ([#164]).

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
