# Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TrixiShallowWater.jl works
with Julia v1.10 and newer. We recommend using the latest stable release of Julia.

## For users
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

## For developers
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