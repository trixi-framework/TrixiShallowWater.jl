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

### For users
TrixiShallowWater.jl and its related tools are registered Julia packages. Hence, you
can install TrixiShallowWater.jl, the numerical solver framework Trixi.jl,
visualization tools [Trixi2Vtk](https://github.com/trixi-framework/Trixi2Vtk.jl), and
[Plots.jl](https://github.com/JuliaPlots/Plots.jl)
as well as relevant time integration sub-packages of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl),
by executing the following commands in the Julia REPL:
```julia
julia> using Pkg

julia> Pkg.add(["TrixiShallowWater", "Trixi", "Trixi2Vtk", "Plots"
                "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqSSPRK"])
```
You can copy and paste all commands to the REPL *including* the leading
`julia>` prompts - they will automatically be stripped away by Julia.
The package [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
and its sub-packages provide time integration schemes used by TrixiShallowWater.jl,
while [Plots.jl](https://github.com/JuliaPlots/Plots.jl) can be used to directly
visualize TrixiShallowWater.jl's results from the REPL.

### For developers
If you plan on editing TrixiShallowWater.jl itself, you can download TrixiShallowWater.jl
locally and use the code from the cloned directory:
```bash
git clone git@github.com:trixi-framework/TrixiShallowWater.jl.git
cd TrixiShallowWater.jl
mkdir run
cd run
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path=".."))' # Install local TrixiShallowWater.jl clone
julia --project=. -e 'using Pkg; Pkg.add(["Trixi", "OrdinaryDiffEq", "Trixi2Vtk", "Plots"])' # Install additional packages
```
Note that the postprocessing tools Trixi2Vtk.jl and Plots.jl are optional and
can be omitted.

If you installed TrixiShallowWater.jl this way, you always have to start Julia with the `--project`
flag set to your `run` directory, e.g.,
```bash
julia --project=.
```
if already inside the `run` directory.
Further details, of how to develop TrixiShallowWater.jl together with a local
clone Trixi.jl, can be found in the Development section of the [documentation](#documentation).

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