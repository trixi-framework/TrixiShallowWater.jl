# Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TrixiShallowWater.jl works
with Julia v1.10 and newer. We recommend using the latest stable release of Julia.

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