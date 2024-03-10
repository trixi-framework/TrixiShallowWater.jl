module TestExamples2DShallowWaterWetDry

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "tree_2d_dgsem")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "TreeMesh2D" begin
#! format: noindent

@testset "Shallow Water Wet/Dry" begin
    @trixi_testset "elixir_shallowwater_ec.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_ec.jl"),
                            l2=[
                                0.991181203601035,
                                0.734130029040644,
                                0.7447696147162621,
                                0.5875351036989047,
                            ],
                            linf=[
                                2.0117744577945413,
                                2.9962317608172127,
                                2.6554999727293653,
                                3.0,
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                0.9130579602987144,
                                1.0602847041965408e-14,
                                1.082225645390032e-14,
                                0.9130579602987147,
                            ],
                            linf=[
                                2.113062037615659,
                                4.6613606802974e-14,
                                5.4225772771633196e-14,
                                2.1130620376156584,
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

    @trixi_testset "elixir_shallowwater_well_balanced_wall.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_wall.jl"),
                            l2=[
                                0.9130579602987144,
                                1.0602847041965408e-14,
                                1.082225645390032e-14,
                                0.9130579602987147,
                            ],
                            linf=[
                                2.113062037615659,
                                4.6613606802974e-14,
                                5.4225772771633196e-14,
                                2.1130620376156584,
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                0.9130579602987147,
                                9.68729463970494e-15,
                                9.694538537436981e-15,
                                0.9130579602987147,
                            ],
                            linf=[
                                2.1130620376156584,
                                2.3875905654916432e-14,
                                2.2492839032269154e-14,
                                2.1130620376156584,
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(flux_lax_friedrichs,
                                                                        hydrostatic_reconstruction_audusse_etal),
                                          flux_nonconservative_audusse_etal),
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with flux_nonconservative_ersing_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                0.9130579602987146,
                                1.0323158914614244e-14,
                                1.0276096319430528e-14,
                                0.9130579602987147,
                            ],
                            linf=[
                                2.11306203761566,
                                4.063916419044386e-14,
                                3.694484044448245e-14,
                                2.1130620376156584,
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_ersing_etal),
                            volume_flux=(flux_wintermeyer_etal,
                                         flux_nonconservative_ersing_etal),
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

    @trixi_testset "elixir_shallowwater_well_balanced_wet_dry.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_wet_dry.jl"),
                            l2=[
                                0.030186039395610056,
                                2.513287752536758e-14,
                                1.3631397744897607e-16,
                                0.10911781485920438,
                            ],
                            linf=[
                                0.49999999999993505,
                                5.5278950497971455e-14,
                                7.462550826772548e-16,
                                2.0,
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

    @trixi_testset "elixir_shallowwater_source_terms.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.001868474306068482,
                                0.01731687445878443,
                                0.017649083171490863,
                                6.274146767717023e-5,
                            ],
                            linf=[
                                0.016962486402209986,
                                0.08768628853889782,
                                0.09038488750767648,
                                0.0001819675955490041,
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

    @trixi_testset "elixir_shallowwater_source_terms_dirichlet.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms_dirichlet.jl"),
                            l2=[
                                0.0018746929418489125,
                                0.017332321628469628,
                                0.01634953679145536,
                                6.274146767717023e-5,
                            ],
                            linf=[
                                0.016262353691956388,
                                0.08726160620859424,
                                0.09043621801418844,
                                0.0001819675955490041,
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_hll" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.0018957692481057034,
                                0.016943229710439864,
                                0.01755623297390675,
                                6.274146767717414e-5,
                            ],
                            linf=[
                                0.015156105797771602,
                                0.07964811135780492,
                                0.0839787097210376,
                                0.0001819675955490041,
                            ],
                            tspan=(0.0, 0.025),
                            surface_flux=(FluxHLL(min_max_speed_naive),
                                          flux_nonconservative_fjordholm_etal))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_nonconservative_ersing_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.002471853426064005,
                                0.05619168608950033,
                                0.11844727575152562,
                                6.274146767730281e-5,
                            ],
                            linf=[
                                0.014332922987500218,
                                0.2141204806174546,
                                0.5392313755637872,
                                0.0001819675955490041,
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_ersing_etal),
                            volume_flux=(flux_wintermeyer_etal,
                                         flux_nonconservative_ersing_etal),
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

    @trixi_testset "elixir_shallowwater_conical_island.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_conical_island.jl"),
                            l2=[
                                0.0459315416430658,
                                0.1644534881916991,
                                0.16445348819169914,
                                0.0011537702354532694,
                            ],
                            linf=[
                                0.21100717610846464,
                                0.9501592344310412,
                                0.9501592344310417,
                                0.021790250683516282,
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
                                0.00025345501281482687,
                                4.4525120338817177e-5,
                                0.00015991819160294247,
                                7.750412064917294e-15,
                            ],
                            linf=[
                                0.004664246019836723,
                                0.0004972780116736669,
                                0.0028735707270457628,
                                6.866729407306593e-14,
                            ],
                            tspan=(0.0, 0.025),
                            basis=LobattoLegendreBasis(3))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_wall.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_wall.jl"),
                            l2=[
                                0.13517233723296504,
                                0.20010876311162215,
                                0.20010876311162223,
                                2.719538414346464e-7,
                            ],
                            linf=[
                                0.5303607982988336,
                                0.5080989745682338,
                                0.5080989745682352,
                                1.1301675764130437e-6,
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
end # SWE

@testset "Two-Layer Shallow Water" begin
    @trixi_testset "elixir_shallowwater_twolayer_convergence.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_twolayer_convergence.jl"),
                            l2=[0.0004016779699408397, 0.005466339651545468,
                                0.006148841330156112,
                                0.0002882339012602492, 0.0030120142442780313,
                                0.002680752838455618,
                                8.873630921431545e-6],
                            linf=[0.002788654460984752, 0.01484602033450666,
                                0.017572229756493973,
                                0.0016010835493927011, 0.009369847995372549,
                                0.008407961775489636,
                                3.361991620143279e-5],
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

    @trixi_testset "elixir_shallowwater_twolayer_well_balanced.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_twolayer_well_balanced.jl"),
                            l2=[3.2935164267930016e-16, 4.6800825611195103e-17,
                                4.843057532147818e-17,
                                0.0030769233188015013, 1.4809161150389857e-16,
                                1.509071695038043e-16,
                                0.0030769233188014935],
                            linf=[2.248201624865942e-15, 2.346382070278936e-16,
                                2.208565017494899e-16,
                                0.026474051138910493, 9.237568031609006e-16,
                                7.520758026187046e-16,
                                0.026474051138910267],
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

    @trixi_testset "elixir_shallowwater_twolayer_well_balanced with flux_lax_friedrichs.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_twolayer_well_balanced.jl"),
                            l2=[2.0525741072929735e-16, 6.000589392730905e-17,
                                6.102759428478984e-17,
                                0.0030769233188014905, 1.8421386173122792e-16,
                                1.8473184927121752e-16,
                                0.0030769233188014935],
                            linf=[7.355227538141662e-16, 2.960836949170518e-16,
                                4.2726562436938764e-16,
                                0.02647405113891016, 1.038795478061861e-15,
                                1.0401789378532516e-15,
                                0.026474051138910267],
                            surface_flux=(flux_lax_friedrichs,
                                          flux_nonconservative_ersing_etal),
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
end # 2LSWE
end # TreeMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
