module TestExamplesTree1D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "tree_1d_dgsem")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "TreeMesh1D" begin
#! format: noindent

@testset "Shallow Water Wet/Dry" begin
    @trixi_testset "elixir_shallowwater_ec.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_ec.jl"),
                            l2=[
                                0.24476140682560343,
                                0.8587309324660326,
                                0.07330427577586297,
                            ],
                            linf=[
                                2.1636963952308372,
                                3.8737770522883115,
                                1.7711213427919539,
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

    @trixi_testset "elixir_shallowwater_ec.jl with initial_condition_weak_blast_wave" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_ec.jl"),
                            l2=[
                                0.39472828074570576,
                                2.0390687947320076,
                                4.1623084150546725e-10,
                            ],
                            linf=[
                                0.7793741954662221,
                                3.2411927977882096,
                                7.419800190922032e-10,
                            ],
                            initial_condition=initial_condition_weak_blast_wave,
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
                                0.10416666834254829,
                                1.4352935256803184e-14,
                                0.10416666834254838,
                            ],
                            linf=[1.9999999999999996, 3.248036646353028e-14, 2.0],
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
                                0.10416666834254835,
                                1.1891029971551825e-14,
                                0.10416666834254838,
                            ],
                            linf=[2.0000000000000018, 2.4019608337954543e-14, 2.0],
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
                                0.10416666834254838,
                                1.6657566141935285e-14,
                                0.10416666834254838,
                            ],
                            linf=[2.0000000000000004, 3.0610625110157164e-14, 2.0],
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
                                0.009657871671690306,
                                6.23878532543316e-14,
                                0.03857583749209947,
                            ],
                            linf=[
                                0.4999999999998892,
                                4.4857500383821557e-13,
                                1.9999999999999811,
                            ],
                            tspan=(0.0, 0.25),
                            # Soften the tolerance as test results vary between different CPUs
                            atol=2000 * eps())
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
                                0.0022363707373868713,
                                0.01576799981934617,
                                4.436491725585346e-5,
                            ],
                            linf=[
                                0.00893601803417754,
                                0.05939797350246456,
                                9.098379777405796e-5,
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
                                0.0022758146627220154,
                                0.015864082886204556,
                                4.436491725585346e-5,
                            ],
                            linf=[
                                0.008457195427364006,
                                0.057201667446161064,
                                9.098379777405796e-5,
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
                                0.005774284062933275,
                                0.017408601639513584,
                                4.43649172561843e-5,
                            ],
                            linf=[
                                0.01639116193303547,
                                0.05102877460799604,
                                9.098379777450205e-5,
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_ersing_etal),
                            volume_flux=(flux_wintermeyer_etal,
                                         flux_nonconservative_ersing_etal),
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
                                0.0022851099219788917,
                                0.01560453773635554,
                                4.43649172558535e-5,
                            ],
                            linf=[
                                0.008934615705174398,
                                0.059403169140869405,
                                9.098379777405796e-5,
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

    @trixi_testset "elixir_shallowwater_source_terms_dirichlet.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms_dirichlet.jl"),
                            l2=[
                                0.0022956052733432287,
                                0.015540053559855601,
                                4.43649172558535e-5,
                            ],
                            linf=[
                                0.008460440313118323,
                                0.05720939349382359,
                                9.098379777405796e-5,
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(FluxHLL(min_max_speed_naive),
                                                                        hydrostatic_reconstruction_audusse_etal),
                                          flux_nonconservative_audusse_etal),
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

    @trixi_testset "elixir_shallowwater_well_balanced_nonperiodic.jl with Dirichlet boundary" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_nonperiodic.jl"),
                            l2=[
                                1.725964362045055e-8,
                                5.0427180314307505e-16,
                                1.7259643530442137e-8,
                            ],
                            linf=[
                                3.844551077492042e-8,
                                3.469453422316143e-15,
                                3.844551077492042e-8,
                            ],
                            tspan=(0.0, 0.25),
                            surface_flux=(FluxHLL(min_max_speed_naive),
                                          flux_nonconservative_fjordholm_etal),)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_well_balanced_nonperiodic.jl with wall boundary" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced_nonperiodic.jl"),
                            l2=[
                                1.7259643614361866e-8,
                                3.5519018243195145e-16,
                                1.7259643530442137e-8,
                            ],
                            linf=[
                                3.844551010878661e-8,
                                9.846474508971374e-16,
                                3.844551077492042e-8,
                            ],
                            tspan=(0.0, 0.25),
                            surface_flux=(FluxHLL(min_max_speed_naive),
                                          flux_nonconservative_fjordholm_etal),
                            boundary_condition=boundary_condition_slip_wall)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_shock_capturing.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_shock_capturing.jl"),
                            l2=[
                                0.07424140641160326,
                                0.2148642632748155,
                                0.0372579849000542,
                            ],
                            linf=[
                                1.1209754279344226,
                                1.3230788645853582,
                                0.8646939843534251,
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

    @trixi_testset "elixir_shallowwater_beach.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_beach.jl"),
                            l2=[
                                0.17979210479598923,
                                1.2377495706611434,
                                6.289818963361573e-8,
                            ],
                            linf=[
                                0.845938394800688,
                                3.3740800777086575,
                                4.4541473087633676e-7,
                            ],
                            tspan=(0.0, 0.05),
                            atol=1e-7) # see https://github.com/trixi-framework/Trixi.jl/issues/1617
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
                                8.965981683033589e-5,
                                1.8565707397810857e-5,
                                4.1043039226164336e-17,
                            ],
                            linf=[
                                0.00041080213807871235,
                                0.00014823261488938177,
                                2.220446049250313e-16,
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
end # SWE

@testset "Two-Layer Shallow Water" begin
    @trixi_testset "elixir_shallowwater_twolayer_convergence.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_twolayer_convergence.jl"),
                            l2=[0.005012009872109003, 0.002091035326731071,
                                0.005049271397924551,
                                0.0024633066562966574, 0.0004744186597732739],
                            linf=[0.0213772149343594, 0.005385752427290447,
                                0.02175023787351349,
                                0.008212004668840978, 0.0008992474511784199],
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
                            l2=[8.949288784402005e-16, 4.0636427176237915e-17,
                                0.001002881985401548,
                                2.133351105037203e-16, 0.0010028819854016578],
                            linf=[2.6229018956769323e-15, 1.878451903240623e-16,
                                0.005119880996670156,
                                8.003199803957679e-16, 0.005119880996670666],
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

    @trixi_testset "elixir_shallowwater_twolayer_dam_break.jl with flux_lax_friedrichs" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_twolayer_dam_break.jl"),
                            l2=[0.1000774903431289, 0.5670692949571057,
                                0.08764242501014498,
                                0.45412307886094555, 0.013638618139749523],
                            linf=[0.586718937495144, 2.1215606128311584,
                                0.5185911311186155,
                                1.820382495072612, 0.5],
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

@testset "Multilayer Shallow Water" begin
    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                0.0012190787287581911,
                                0.0007900811065315603,
                                0.00024265888138436312,
                                0.00044695521604641624,
                                0.0009421189161161149,
                                0.0001681269882698457,
                                0.000139126377287091,
                            ],
                            linf=[
                                0.0035625489595796367,
                                0.0025099269565310167,
                                0.000677952137713711,
                                0.0019817824651813254,
                                0.002749046992753801,
                                0.0006393958093637853,
                                0.00021874455861881081,
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

    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl with flux_lax_friedrichs" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                0.0007413423106723253,
                                0.000540576469967569,
                                8.772494905937258e-5,
                                0.0006969877681271723,
                                0.0004924792705658947,
                                0.00011816401639126746,
                                0.000139126377287091,
                            ],
                            linf=[
                                0.0026570823889522366,
                                0.0028102271323960926,
                                0.0002285569562520129,
                                0.003310623577850169,
                                0.0024909798619617285,
                                0.0005211057995907487,
                                0.00021874455861881081,
                            ],
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

    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl with hydrostatic_reconstruction_ersing_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                0.0007413415478415575,
                                0.0005405759470722189,
                                8.772503316220759e-5,
                                0.0006969870846204186,
                                0.0004924774643042689,
                                0.00011816402423681591,
                                0.000139126377287091,
                            ],
                            linf=[
                                0.0026571096352370205,
                                0.0028102212809490434,
                                0.00022855702706781056,
                                0.0033106637858648646,
                                0.0024909780761511735,
                                0.0005211084155695156,
                                0.00021874455861881081,
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                                            DissipationLocalLaxFriedrichs()),
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
                                1.2015338985485869e-17,
                                1.4450281975802518e-17,
                                0.001002881985401688,
                                5.6968770683622844e-18,
                                5.8630384557660415e-18,
                                1.938143907661441e-17,
                                0.0010028819854016856,
                            ],
                            linf=[
                                8.326672684688674e-17,
                                1.3877787807814457e-16,
                                0.005119880996670767,
                                3.2043414465411084e-17,
                                3.288278489632501e-17,
                                1.6336799907854548e-16,
                                0.005119880996670708,
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced_wet_dry.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced_wet_dry.jl"),
                            l2=[
                                0.025515633605684252,
                                0.020321228314358498,
                                2.078418885421747e-16,
                                1.5034457079135132e-17,
                                0.04746160729122716,
                            ],
                            linf=[
                                0.9999999999999989,
                                1.0,
                                1.8096275952411747e-15,
                                4.245490280714591e-16,
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced_wet_dry.jl with perturbation" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced_wet_dry.jl"),
                            l2=[
                                0.031114402256402142,
                                0.02924759823866292,
                                0.0021533928763101534,
                                0.01412031518684538,
                                0.04746160729122716,
                            ],
                            linf=[
                                1.2499999999999991,
                                0.998672198687742,
                                0.04297148454930247,
                                0.24675612377725484,
                                2.0,
                            ],
                            tspan=(0.0, 0.25),
                            perturbation=true)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_multilayer_dam_break_ec.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break_ec.jl"),
                            l2=[
                                0.03797459518473401,
                                0.038418485935199,
                                0.051223819840353284,
                                0.29104239802837284,
                                0.2897819870634873,
                                0.21991790700220484,
                                0.013638618139749504,
                            ],
                            linf=[
                                0.14811281751508787,
                                0.14294051443354894,
                                0.7862151912522912,
                                1.079364356014889,
                                1.0144419583593887,
                                0.992467173339586,
                                0.4999999999999993,
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break_es.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break_es.jl"),
                            l2=[
                                0.037157941241555185,
                                0.03944955271506525,
                                0.0473661084974011,
                                0.28686995762998163,
                                0.2869427613009073,
                                0.21860712413448047,
                                0.013638618139749504,
                            ],
                            linf=[
                                0.18060907942891546,
                                0.2221666684315402,
                                0.44639752069288585,
                                1.146783259827299,
                                1.140481663033156,
                                0.9997998219394365,
                                0.4999999999999993,
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break_es_dry.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break_es_dry.jl"),
                            l2=[
                                0.0539713484191238,
                                0.05378401084239003,
                                0.05343571027541134,
                                0.05294334126567629,
                                0.14409266517044308,
                                0.2258690603183475,
                                0.22395024283599815,
                                0.22041584295217526,
                                0.21549254848916716,
                                0.6274449804372626,
                                0.013638618139749523,
                            ],
                            linf=[
                                0.27418120077237185,
                                0.2738885780879326,
                                0.2733530772398547,
                                0.2726148157344642,
                                0.8243308057901404,
                                0.6976533131655047,
                                0.6927619903064818,
                                0.6836813560099443,
                                0.6708713913679563,
                                3.0233275007119254,
                                0.5,
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

    @trixi_testset "elixir_shallowwater_multilayer_beach.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_beach.jl"),
                            l2=[
                                0.5837610281036327,
                                3.7410178422784273,
                                6.289710784508112e-8,
                            ],
                            linf=[
                                0.9917260088977286,
                                7.0962810662023035,
                                4.452604027704865e-7,
                            ],
                            tspan=(0.0, 0.25),
                            atol=1e-5)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_multilayer_parabolic_bowl.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_parabolic_bowl.jl"),
                            l2=[
                                0.002063065256405533,
                                0.00048255459387538827,
                                4.0617503413126664e-17,
                            ],
                            linf=[
                                0.0036250458403655466,
                                0.0008665322585535845,
                                1.6653345369377348e-16,
                            ],
                            tspan=(0.0, 0.25),
                            atol=1e-9,
                            coverage_override=(maxiters = 15,))
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
end # TreeMesh1D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
