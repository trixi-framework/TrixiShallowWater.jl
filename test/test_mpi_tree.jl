module TestExamplesMPITreeMesh

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

const EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "tree_2d_dgsem")

@testset "TreeMesh MPI" begin
#! format: noindent

@testset "Shallow Water Wet/Dry" begin
    @trixi_testset "elixir_shallowwater_ec.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_ec.jl"),
                            l2=[
                                0.9911802019934329,
                                0.7340106828033273,
                                0.7446338002084801,
                                0.5875351036989047
                            ],
                            linf=[
                                2.0120253138457564,
                                2.991158989293406,
                                2.6557412817714035,
                                3.0
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
                                0.9130579602987147
                            ],
                            linf=[
                                2.113062037615659,
                                4.6613606802974e-14,
                                5.4225772771633196e-14,
                                2.1130620376156584
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

    @trixi_testset "elixir_shallowwater_source_terms_dirichlet.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms_dirichlet.jl"),
                            l2=[
                                0.0018596727473552813,
                                0.017306217777629147,
                                0.016367646997420396,
                                6.274146767723934e-5
                            ],
                            linf=[
                                0.016548007102923368,
                                0.08726160568822783,
                                0.09043621622245013,
                                0.0001819675955490041
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

    @trixi_testset "elixir_shallowwater_conical_island.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_conical_island.jl"),
                            l2=[
                                0.045928568956367245,
                                0.16446498697148945,
                                0.16446498697148945,
                                0.0011537702354532694
                            ],
                            linf=[
                                0.21098104635388315,
                                0.9501826412445212,
                                0.9501826412445218,
                                0.021790250683516282
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
                                0.0002534446645366419,
                                4.452398574250411e-5,
                                0.0001599158029172979,
                                1.0761492186651631e-16
                            ],
                            linf=[
                                0.00466424300648733,
                                0.0004972785186291824,
                                0.0028735682830073137,
                                6.661338147750939e-16
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

    @trixi_testset "elixir_shallowwater_inflow_outflow.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_inflow_outflow.jl"),
                            l2=[
                                0.164716617086721,
                                0.5257126140039803,
                                0.5257126140039803,
                                0.0
                            ],
                            linf=[
                                0.5595760580954796,
                                1.3874204364467229,
                                1.3874204364467246,
                                0.0
                            ],
                            # Increase iterations for coverage testing to trigger both inflow and outflow conditions
                            coverage_override=(maxiters = 65, tspan = (0.0, 2.5)))
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
end # 2LSWE

@testset "Multilayer Shallow Water" begin
    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                0.001578480912165503,
                                0.0010725791806397446,
                                0.00036117546296274316,
                                0.0010107766709183063,
                                0.0009612673392275463,
                                0.00026390576962071974,
                                0.001140711809558893,
                                0.001197469124790338,
                                0.00031384421384170294,
                                0.00019675440964325044
                            ],
                            linf=[
                                0.006611273262461914,
                                0.0047987596637860674,
                                0.001476163166610922,
                                0.004378481797734368,
                                0.004059620398096986,
                                0.00107253722812789,
                                0.004650692998510841,
                                0.005014864957408438,
                                0.0013175223929244861,
                                0.0004374891172380657
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break.jl"),
                            l2=[
                                0.006867196239306527,
                                0.007212095601190101,
                                0.014309127811972761,
                                0.005801387118671292,
                                0.0057555558162007536,
                                0.006951246408270715,
                                6.084817535430214e-5,
                                6.068733435363927e-5,
                                7.807791370377477e-5,
                                0.003857583749542185
                            ],
                            linf=[
                                0.028792046975677804,
                                0.036503344705121676,
                                0.19267105917517297,
                                0.028685238833612656,
                                0.02761561086221269,
                                0.029892140357308587,
                                0.0004158435510054702,
                                0.0004119057098728111,
                                0.0007709896609887244,
                                0.10000011323773067
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
end # MLSWE
end # TreeMesh MPI

end # module
