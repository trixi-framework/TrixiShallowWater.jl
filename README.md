# TrixiShallowWater.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://trixi-framework.github.io/TrixiShallowWater.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://trixi-framework.github.io/TrixiShallowWater.jl/dev/)
[![Build Status](https://github.com/trixi-framework/TrixiShallowWater.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/trixi-framework/TrixiShallowWater.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/trixi-framework/TrixiShallowWater.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/trixi-framework/TrixiShallowWater.jl)
[![Coverage](https://coveralls.io/repos/github/trixi-framework/TrixiShallowWater.jl/badge.svg?branch=main)](https://coveralls.io/github/trixi-framework/TrixiShallowWater.jl?branch=main)

**TrixiShallowWater.jl** is a numerical simulation package focused on solving shallow water equations 
with the discontinuous Galerkin method and written in Julia. The package builds on the numerical solver [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) 
and provides several specialized models and features specific for shallow water applications.

## Examples

## Installation
TrixiShallowWater.jl is **not** a registered Julia package, and therefore needs to be downloaded manually and then run from with the cloned directory:
```bash
git clone https://github.com/trixi-framework/TrixiShallowWater.jl.git
julia --project=@.
```
In addition TrixiShallowWater.jl requires the numerical solver framework [Trixi.jl](https://github.com/trixi-framework/Trixi.jl), [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) for time integration, and [Plots.jl](https://github.com/JuliaPlots/Plots.jl) for visualization, which can be installed by executing the following in the Julia REPL:
```julia
julia> using Pkg

julia> Pkg.add(["Trixi", "Trixi2Vtk", "OrdinaryDiffEq", "Plots"])
```

## License and contributing
TrixiShallowWater.jl is licensed under the MIT license (see [LICENSE.md](LICENSE.md)). Since Trixi.jl is
an open-source project, we are very happy to accept contributions from the
community. To get in touch with the developers,
[join us on Slack](https://join.slack.com/t/trixi-framework/shared_invite/zt-sgkc6ppw-6OXJqZAD5SPjBYqLd8MU~g)
or [create an issue](https://github.com/trixi-framework/TrixiShallowWater.jl/issues/new).


**Note: This repository is still in its alpha stage and anything might change at
any time and without warning, including the deletion of this repository
itself.**
