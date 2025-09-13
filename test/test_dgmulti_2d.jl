module TestExamplesDGMulti2D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "dgmulti_2d")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "DGMulti 2D" begin
#! format: noindent

@trixi_testset "elixir_shallowwater_source_terms.jl (Quad, SBP)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_source_terms.jl"),
                        cells_per_dimension=8, element_type=Quad(),
                        approximation_type=SBP(),
                        l2=[
                            0.0020316463892983217,
                            0.02366902012965938,
                            0.03446194535725363,
                            1.921676942941478e-15
                        ],
                        linf=[
                            0.010384996665098178,
                            0.08750632767286826,
                            0.12088391569555768,
                            9.325873406851315e-15
                        ])
    @test_allocations(Trixi.rhs!, semi, sol, 1000)
end

@trixi_testset "elixir_shallowwater_source_terms.jl (Tri, SBP)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_source_terms.jl"),
                        cells_per_dimension=8, element_type=Tri(),
                        approximation_type=SBP(),
                        l2=[
                            0.004180679992535108,
                            0.07026193567927695,
                            0.11815151184746633,
                            2.3786840926019625e-15
                        ],
                        linf=[
                            0.020760033097378283,
                            0.29169608872805686,
                            0.567418412384793,
                            1.1102230246251565e-14
                        ])
    @test_allocations(Trixi.rhs!, semi, sol, 1000)
end

@trixi_testset "elixir_shallowwater_source_terms.jl (Tri, Polynomial)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_source_terms.jl"),
                        cells_per_dimension=8, element_type=Tri(),
                        approximation_type=Polynomial(),
                        # The last l2, linf error are the L2 projection error in approximating `b`, so they are not
                        # zero for general non-collocated quadrature rules (e.g., for `element_type=Tri()`, `polydeg > 2`).
                        l2=[
                            0.0008309358577296097,
                            0.015224511207450263,
                            0.016033971785878454,
                            1.282024730815488e-5
                        ],
                        linf=[
                            0.0018880416154898327,
                            0.05466845626696504,
                            0.06345896594568323,
                            3.398993309877696e-5
                        ])
    @test_allocations(Trixi.rhs!, semi, sol, 1000)
end

@trixi_testset "elixir_shallowwater_source_terms.jl (Quad, Polynomial)" begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_source_terms.jl"),
                        cells_per_dimension=8, element_type=Quad(),
                        approximation_type=Polynomial(),
                        # The last l2, linf error are the L2 projection error in approximating `b`. However, this is zero
                        # for `Quad()` elements with `Polynomial()` approximations because the quadrature rule defaults to
                        # a `(polydeg + 1)`-point Gauss quadrature rule in each coordinate (in general, StartUpDG.jl defaults
                        # to the quadrature rule with the fewest number of points which exactly integrates the mass matrix).
                        l2=[
                            7.460473151203597e-5,
                            0.0036855901000765463,
                            0.003910160802530521,
                            6.743418333559633e-15
                        ],
                        linf=[
                            0.0002599957400737374,
                            0.007223608258381642,
                            0.010364657535841815,
                            2.042810365310288e-14
                        ])
    @test_allocations(Trixi.rhs!, semi, sol, 1000)
end
end # DGMulti2D

end # module
