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
                                0.244729018751225,
                                0.8583565222389505,
                                0.07330427577586297,
                            ],
                            linf=[
                                2.1635021283528504,
                                3.8717508164234453,
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
                                0.39464782107209717,
                                2.03880864210846,
                                4.1623084150546725e-10,
                            ],
                            linf=[
                                0.778905801278281,
                                3.2409883402608273,
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
                                0.00965787167169024,
                                5.345454081916856e-14,
                                0.03857583749209928,
                            ],
                            linf=[
                                0.4999999999998892,
                                2.2447689894899726e-13,
                                1.9999999999999714,
                            ],
                            tspan=(0.0, 0.25),
                            # Soften the tolerance as test results vary between different CPUs
                            atol=1000 * eps())
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
                                0.005012009872110106,
                                0.005049271397926405,
                                0.0020910353267268073,
                                0.002463306656296269,
                                0.0004744186597731183,
                            ],
                            linf=[
                                0.021377214934391375,
                                0.021750237873536582,
                                0.005385752427245816,
                                0.008212004668841977,
                                0.0008992474511775317,
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
end # MLSWE
end # TreeMesh1D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
