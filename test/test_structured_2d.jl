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
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
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
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
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
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_conical_island.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_conical_island.jl"),
                            l2=[
                                0.04593154164306353,
                                0.1644534881916908,
                                0.16445348819169076,
                                0.0011537702354532122
                            ],
                            linf=[
                                0.21100717610846442,
                                0.9501592344310412,
                                0.950159234431041,
                                0.021790250683516296
                            ],
                            tspan=(0.0, 0.025))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_parabolic_bowl.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_parabolic_bowl.jl"),
                            l2=[
                                0.00015285369980313484,
                                1.9536806395943226e-5,
                                9.936906607758672e-5,
                                5.0686313334616055e-15
                            ],
                            linf=[
                                0.003316119030459211,
                                0.0005075409427972817,
                                0.001986721761060583,
                                4.701794509287538e-14
                            ],
                            tspan=(0.0, 0.025), cells_per_dimension=(40, 40))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end
end # SWE
end # StructuredMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
