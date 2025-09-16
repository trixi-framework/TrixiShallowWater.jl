module TestExamplesT8codeMesh2D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "t8code_2d_dgsem")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "T8codeMesh2D" begin
#! format: noindent

@trixi_testset "elixir_shallowwater_source_terms.jl" begin
    # This test is identical to the one in `test_p4est_2d.jl`.
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_source_terms.jl"),
                        l2=[
                            9.168126407325352e-5,
                            0.0009795410115453788,
                            0.002546408320320785,
                            3.941189812642317e-6
                        ],
                        linf=[
                            0.0009903782521019089,
                            0.0059752684687262025,
                            0.010941106525454103,
                            1.2129488214718265e-5
                        ],
                        tspan=(0.0, 0.1))
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(Trixi.rhs!, semi, sol, 1000)
end
end # T8codeMesh2D
end # module
