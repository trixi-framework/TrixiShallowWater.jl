module TestExamplesStructuredMesh2D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "structured_2d_dgsem")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "StructuredMesh2D" begin
#! format: noindent

@testset "Shallow Water Wet/Dry" begin
    @trixi_testset "elixir_shallowwater_source_terms.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.0017286908591070864,
                                0.025585037307655684,
                                0.028374244567802766,
                                6.274146767730866e-5
                            ],
                            linf=[
                                0.012973752001194772,
                                0.10829375385832263,
                                0.15832858475438094,
                                0.00018196759554722775
                            ],
                            tspan=(0.0, 0.05))
        @test_allocations(Trixi.rhs!, semi, sol,  1000)
    end

    @trixi_testset "elixir_shallowwater_well_balanced.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                0.7920927046419308,
                                9.92129670988898e-15,
                                1.0118635033124588e-14,
                                0.7920927046419308
                            ],
                            linf=[
                                2.408429868800133,
                                5.5835419986809516e-14,
                                5.448874313931364e-14,
                                2.4084298688001335
                            ],
                            tspan=(0.0, 0.25))
        @test_allocations(Trixi.rhs!, semi, sol,  1000)
    end

    @trixi_testset "elixir_shallowwater_well_balanced_wet_dry.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_wet_dry.jl"),
                            l2=[
                                0.019731646454942086,
                                1.0694532773278277e-14,
                                1.1969913383405568e-14,
                                0.0771517260037954
                            ],
                            linf=[
                                0.4999999999998892,
                                6.067153702623552e-14,
                                4.4849667259339357e-14,
                                1.9999999999999993
                            ],
                            tspan=(0.0, 0.25),
                            atol=1e-12)
        @test_allocations(Trixi.rhs!, semi, sol,  1000)
    end

    @trixi_testset "elixir_shallowwater_conical_island.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_conical_island.jl"),
                            l2=[
                                0.04592856895636503,
                                0.16446498697148132,
                                0.16446498697148126,
                                0.0011537702354532122
                            ],
                            linf=[
                                0.21098104635388404,
                                0.950182641244522,
                                0.950182641244521,
                                0.021790250683516296
                            ],
                            tspan=(0.0, 0.025))
        @test_allocations(Trixi.rhs!, semi, sol,  1000)
    end

    @trixi_testset "elixir_shallowwater_parabolic_bowl.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_parabolic_bowl.jl"),
                            l2=[
                                0.00015285044801809777,
                                1.953670773024699e-5,
                                9.937604367810603e-5,
                                9.676656822027896e-17
                            ],
                            linf=[
                                0.0033161186020839208,
                                0.0005075478178365348,
                                0.001986721591574139,
                                7.771561172376096e-16
                            ],
                            tspan=(0.0, 0.025), cells_per_dimension=(40, 40))
        @test_allocations(Trixi.rhs!, semi, sol,  1000)
    end
end # SWE
end # StructuredMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
