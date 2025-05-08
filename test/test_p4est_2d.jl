module TestExamplesP4estMesh2D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "p4est_2d_dgsem")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "P4estMesh2D" begin
#! format: noindent

@testset "Shallow Water Wet/Dry" begin
    @trixi_testset "elixir_shallowwater_source_terms.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
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
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_well_balanced_nonconforming.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_nonconforming.jl"),
                            l2=[
                                0.2018723974651268,
                                4.798334932564209e-14,
                                5.265746106509223e-14,
                                0.20187239746512692
                            ],
                            linf=[
                                0.41600528180178653,
                                1.4686699343325025e-12,
                                2.6443321644596183e-12,
                                0.41600528180178625
                            ],
                            tspan=(0.0, 0.25),
                            atol=1e-10)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_perturbation_amr.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_perturbation_amr.jl"),
                            l2=[
                                0.02263230105470324,
                                0.09090425233020173,
                                0.09124622065757255,
                                0.0011045848311422332
                            ],
                            linf=[
                                0.3118823726810007,
                                0.7855402508435719,
                                0.7401368273982774,
                                0.011669083581857587
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

    @trixi_testset "elixir_shallowwater_well_balanced_wet_dry_nonconforming.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_wet_dry_nonconforming.jl"),
                            l2=[
                                0.17389058166483437,
                                1.6059906261139574e-14,
                                1.7488462093286073e-14,
                                0.18415049792609273
                            ],
                            linf=[
                                0.4160052818017864,
                                5.293485363361597e-13,
                                9.214443365763228e-13,
                                0.41600528180178625
                            ],
                            tspan=(0.0, 0.25),
                            atol=1e-10)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    # Note, these values may change as the functionality of well-balanced mortars
    # with AMR and wet/dry are further developed according to the issue
    # https://github.com/trixi-framework/TrixiShallowWater.jl/issues/77
    @trixi_testset "elixir_shallowwater_perturbation_wet_dry_amr.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_perturbation_wet_dry_amr.jl"),
                            l2=[
                                0.399907298565006,
                                0.3899824006669886,
                                0.44130722321781896,
                                0.4564350570616199
                            ],
                            linf=[
                                1.3905078647694498,
                                3.3601802680830555,
                                3.576650520299362,
                                0.7495177590247986
                            ],
                            tspan=(0.0, 0.025),
                            coverage_override=(maxiters = 5,))
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
end # P4estMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
