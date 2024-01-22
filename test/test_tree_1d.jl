module TestExamplesTree1D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "tree_1d_dgsem")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "TreeMesh1D" begin
#! format: noindent

# Run basic tests
@testset "Examples 1D" begin
    # Shallow water
    include("test_tree_1d_shallowwater_wet_dry.jl")
end

# Clean up afterwards: delete Trixi.jl output directory
@test_nowarn rm(outdir, recursive = true)
end # TreeMesh1D

end # module
