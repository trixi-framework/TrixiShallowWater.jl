# TrixiShallowWater.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://trixi-framework.github.io/TrixiShallowWater.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://trixi-framework.github.io/TrixiShallowWater.jl/dev/)
[![Slack](https://img.shields.io/badge/chat-slack-e01e5a)](https://join.slack.com/t/trixi-framework/shared_invite/zt-sgkc6ppw-6OXJqZAD5SPjBYqLd8MU~g)
[![Build Status](https://github.com/trixi-framework/TrixiShallowWater.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/trixi-framework/TrixiShallowWater.jl/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/trixi-framework/TrixiShallowWater.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/trixi-framework/TrixiShallowWater.jl)
[![Coverage](https://coveralls.io/repos/github/trixi-framework/TrixiShallowWater.jl/badge.svg?branch=main)](https://coveralls.io/github/trixi-framework/TrixiShallowWater.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15206520.svg)](https://doi.org/10.5281/zenodo.15206520)

<p align="center">
  <img width="300px" src="https://trixi-framework.github.io/assets/logo_sw.png">
</p>

**TrixiShallowWater.jl** is a numerical simulation package focused on solving shallow water equations
with the discontinuous Galerkin method and written in Julia. The package builds on the numerical
simulation framework for conservation laws [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
and provides several specialized models and features specific for shallow water applications.
Below is a short summary of the available features:

* 1D and 2D simulations on [line/quad meshes](https://trixi-framework.github.io/Trixi.jl/stable/overview/#Semidiscretizations)
  * Cartesian and curvilinear meshes
  * Conforming and non-conforming meshes
  * Hierarchical quadtree meshes with adaptive mesh refinement
* High-order accuracy in space and time
  * Entropy-stable discontinuous Galerkin methods based on flux differencing
  * Entropy-stable shock capturing
  * Positivity-preserving limiting
  * Compatible with the [SciML ecosystem for ordinary differential equations](https://diffeq.sciml.ai/latest/)
  * CFL-based and error-based time step control
* Shallow water capabilities
  * Wetting and drying
  * Multi-layer flows
  * Sediment transport via an Exner model

## Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TrixiShallowWater.jl works
with Julia v1.10 and newer. We recommend using the latest stable release of Julia.

[comment]: <> (We can update this with a "for users" and "for developers" section once the package is registered)

TrixiShallowWater.jl is **not** a registered Julia package, and therefore needs to be downloaded manually and then run from within the cloned directory:
```bash
git clone https://github.com/trixi-framework/TrixiShallowWater.jl.git
julia --project=@.
```
In addition TrixiShallowWater.jl requires the numerical solver framework [Trixi.jl](https://github.com/trixi-framework/Trixi.jl), relevant sub-packages of [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) for time integration, and [Plots.jl](https://github.com/JuliaPlots/Plots.jl) for visualization, which can be installed by executing the following in the Julia REPL:
```julia
julia> using Pkg

julia> Pkg.add(["Trixi", "Trixi2Vtk", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqSSPRK", "Plots"])
```

## Authors
TrixiShallowWater.jl is maintained by the
[Trixi authors](https://github.com/trixi-framework/Trixi.jl/blob/main/AUTHORS.md).
Its principal developers are [Andrew Winters](https://liu.se/en/employee/andwi94)
(Linköping University, Sweden)
and [Patrick Ersing](https://liu.se/en/employee/pater53)
(Linköping University, Sweden).
The full list of contributors can be found in [AUTHORS.md](AUTHORS.md).

## License and contributing
TrixiShallowWater.jl is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
Since TrixiShallowWater.jl is an open-source project, we are very happy to accept contributions from the
community. To get in touch with the developers,
[join us on Slack](https://join.slack.com/t/trixi-framework/shared_invite/zt-sgkc6ppw-6OXJqZAD5SPjBYqLd8MU~g)
or [create an issue](https://github.com/trixi-framework/TrixiShallowWater.jl/issues/new).

## Acknowledgments
<p align="center" style="font-size:0;"><!--
  SRC      --><img align="middle" src="https://github.com/trixi-framework/Trixi.jl/assets/3637659/48f9da06-6f7a-4586-b23e-739bee3901c0" height="120"><!--
  -->
</p>

This project has benefited from funding from [Vetenskapsrådet](https://www.vr.se)
(VR, Swedish Research Council), Sweden
through the VR Starting Grant "Shallow water flows including sediment transport and morphodynamics",
VR grant agreement 2020-03642 VR.