module TestExamples2DShallowWater

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

    @trixi_testset "elixir_shallowwater_well_balanced_wall.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_wall.jl"),
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                0.9130579602987147,
                                9.68729463970494e-15,
                                9.694538537436981e-15,
                                0.9130579602987147
                            ],
                            linf=[
                                2.1130620376156584,
                                2.3875905654916432e-14,
                                2.2492839032269154e-14,
                                2.1130620376156584
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(FluxLaxFriedrichs(max_abs_speed_naive),
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with flux_nonconservative_wintermeyer_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                0.9130579602987146,
                                1.0323158914614244e-14,
                                1.0276096319430528e-14,
                                0.9130579602987147
                            ],
                            linf=[
                                2.11306203761566,
                                4.063916419044386e-14,
                                3.694484044448245e-14,
                                2.1130620376156584
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_wintermeyer_etal),
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
                                0.10911781485920438
                            ],
                            linf=[
                                0.49999999999993505,
                                5.5278950497971455e-14,
                                7.462550826772548e-16,
                                2.0
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
                                6.274146767717023e-5
                            ],
                            linf=[
                                0.016962486402209986,
                                0.08768628853889782,
                                0.09038488750767648,
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_hll" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.0018957692481057034,
                                0.016943229710439864,
                                0.01755623297390675,
                                6.274146767717414e-5
                            ],
                            linf=[
                                0.015156105797771602,
                                0.07964811135780492,
                                0.0839787097210376,
                                0.0001819675955490041
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with FluxHLL(min_max_speed_naive)" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.0018957692481057034,
                                0.016943229710439864,
                                0.01755623297390675,
                                6.274146767717414e-5
                            ],
                            linf=[
                                0.015156105797771602,
                                0.07964811135780492,
                                0.0839787097210376,
                                0.0001819675955490041
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_nonconservative_wintermeyer_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.002471853426064005,
                                0.05619168608950033,
                                0.11844727575152562,
                                6.274146767730281e-5
                            ],
                            linf=[
                                0.014332922987500218,
                                0.2141204806174546,
                                0.5392313755637872,
                                0.0001819675955490041
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_wintermeyer_etal),
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with dissipation_roe" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.00032946467230644263,
                                0.017644710574401774,
                                0.04741010918138461,
                                6.274146767723582e-5
                            ],
                            linf=[
                                0.0021016896757704018,
                                0.06358320986874899,
                                0.15285461591374805,
                                0.0001819675955490041
                            ],
                            surface_flux=(FluxPlusDissipation(flux_central,
                                                              dissipation_roe),
                                          flux_nonconservative_fjordholm_etal),
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

    @trixi_testset "elixir_shallowwater_wall.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_wall.jl"),
                            l2=[
                                0.1351723240085936,
                                0.20010881416550014,
                                0.2001088141654999,
                                2.719538414346464e-7
                            ],
                            linf=[
                                0.5303608302490757,
                                0.5080987791967457,
                                0.5080987791967506,
                                1.1301675764130437e-6
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

    @trixi_testset "elixir_shallowwater_inflow_outflow_reverse.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_inflow_outflow.jl"),
                            l2=[
                                0.16471661708672095,
                                0.52571261400398,
                                0.5257126140039801,
                                0.0
                            ],
                            linf=[
                                0.5595760580954816,
                                1.3874204364467226,
                                1.3874204364467244,
                                0.0
                            ],
                            v1=0.1, v2=0.1,
                            boundary_condition_inflow=BoundaryConditionMomentum(t -> 0.1 -
                                                                                     0.05 *
                                                                                     t,
                                                                                t -> 0.1 -
                                                                                     0.05 *
                                                                                     t,
                                                                                equations),
                            boundary_conditions=(x_neg = boundary_condition_outflow,
                                                 x_pos = boundary_condition_inflow,
                                                 y_neg = boundary_condition_outflow,
                                                 y_pos = boundary_condition_inflow),
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

    @trixi_testset "elixir_shallowwater_twolayer_well_balanced with FluxLaxFriedrichs(max_abs_speed_naive).jl" begin
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
                            surface_flux=(FluxLaxFriedrichs(max_abs_speed_naive),
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

    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl with FluxLaxFriedrichs(max_abs_speed_naive)" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                0.0009657284251663235,
                                0.0007574316192008629,
                                0.00014884847789419973,
                                0.0008365221427534766,
                                0.0005829573732731771,
                                0.00013759996082561183,
                                0.0009795167874289885,
                                0.0007208939038902295,
                                0.00015863343571957187,
                                0.00019675440964325044
                            ],
                            linf=[
                                0.005343578840408814,
                                0.005546044302341735,
                                0.0005110705217063471,
                                0.004897047461771331,
                                0.004047517471484768,
                                0.0005991894625433924,
                                0.006191481703377466,
                                0.005110806058533313,
                                0.000775668423696807,
                                0.0004374891172380657
                            ],
                            surface_flux=(FluxLaxFriedrichs(max_abs_speed_naive),
                                          flux_nonconservative_ersing_etal),
                            tspan=(0.0, 0.25),
                            atol=1e-11)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                0.0015784809121655386,
                                0.0010725791806396989,
                                0.000361175462962745,
                                0.001010776670918376,
                                0.0009612673392275332,
                                0.0002639057696207499,
                                0.0011407118095590378,
                                0.0011974691247903164,
                                0.00031384421384170354,
                                0.00019675440964325044
                            ],
                            linf=[
                                0.006611273262456141,
                                0.004798759663784402,
                                0.0014761631666105335,
                                0.004378481797736367,
                                0.004059620398099983,
                                0.001072537228128334,
                                0.004650692998510397,
                                0.00501486495741188,
                                0.0013175223929257074,
                                0.0004374891172380657
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(flux_ersing_etal,
                                                                        hydrostatic_reconstruction_ersing_etal),
                                          FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                                                        hydrostatic_reconstruction_ersing_etal)),
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced.jl"),
                            l2=[
                                8.476558640011521e-18,
                                8.37932362312502e-18,
                                1.6048527438048627e-17,
                                0.003076923318801488,
                                4.686259874366214e-18,
                                4.2109008869298085e-18,
                                3.738959845788582e-18,
                                2.202609051659983e-17,
                                4.6862598743662146e-18,
                                4.2109008869298085e-18,
                                3.738959845788581e-18,
                                2.2026090516599833e-17,
                                0.0030769233188014935
                            ],
                            linf=[
                                1.8041124150158794e-16,
                                7.632783294297951e-17,
                                2.3592239273284576e-16,
                                0.026474051138910215,
                                6.685363189156086e-17,
                                5.2914444984119546e-17,
                                6.403718148497213e-17,
                                2.653446487129347e-16,
                                6.685363189156086e-17,
                                5.2914444984119546e-17,
                                6.403718148497213e-17,
                                2.653446487129347e-16,
                                0.026474051138910267
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced.jl with FluxLaxFriedrichs(max_abs_speed_naive)" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced.jl"),
                            l2=[
                                5.2242434650072826e-18,
                                7.520718359606061e-18,
                                6.0112387963921125e-18,
                                0.0030769233188014883,
                                5.827683385814391e-18,
                                6.3826870710212345e-18,
                                5.614916422687807e-18,
                                3.6331769798446284e-17,
                                5.8276833858143906e-18,
                                6.382687071021235e-18,
                                5.614916422687807e-18,
                                3.633176979844629e-17,
                                0.0030769233188014935
                            ],
                            linf=[
                                3.469446951953614e-17,
                                4.85722573273506e-17,
                                5.551115123125783e-17,
                                0.026474051138910215,
                                3.440924806544217e-17,
                                4.517518549935627e-17,
                                4.188211194330907e-17,
                                2.034006113871492e-16,
                                3.440924806544217e-17,
                                4.517518549935627e-17,
                                4.188211194330907e-17,
                                2.034006113871492e-16,
                                0.026474051138910267
                            ],
                            surface_flux=(FluxLaxFriedrichs(max_abs_speed_naive),
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break with FluxLaxFriedrichs(max_abs_speed_naive).jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break.jl"),
                            l2=[
                                0.007162130568658371,
                                0.007507487998680027,
                                0.01500636533827742,
                                0.005593279150854283,
                                0.005603720137202845,
                                0.0068519457928639385,
                                5.748086634734491e-5,
                                5.7787474070077486e-5,
                                7.492163249697633e-5,
                                0.003857583749542185
                            ],
                            linf=[
                                0.037711432402825124,
                                0.0427565001934887,
                                0.17275386567546575,
                                0.025348821041845597,
                                0.02434895938442919,
                                0.03396954814003092,
                                0.0003026537474088657,
                                0.00029947595638190504,
                                0.00037146307629417057,
                                0.10000011323773067
                            ],
                            surface_flux=(FluxLaxFriedrichs(max_abs_speed_naive),
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break_dry.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break_dry.jl"),
                            l2=[
                                0.030684387782428005,
                                0.03094529104416787,
                                0.09023105191044616,
                                0.021006479937144513,
                                0.02094556081555727,
                                0.06145762035878689,
                                0.0010401307741441598,
                                0.0010346793974473935,
                                0.002952821531601029,
                                0.005470808030402103
                            ],
                            linf=[
                                0.10807087634534615,
                                0.11202631201715663,
                                0.33977650047095415,
                                0.05694587797245012,
                                0.056598515593249694,
                                0.17171444846932465,
                                0.00827507857831736,
                                0.008095427727762589,
                                0.026626783423061535,
                                0.1016120899921184
                            ],
                            tspan=(0.0, 0.25),
                            # Increase iterations for coverage testing to trigger the
                            # positivity limiter
                            coverage_override=(maxiters = 130, tspan = (0.0, 1.5)))
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
end # TreeMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
