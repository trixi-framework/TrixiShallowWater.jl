module TestExamplesUnstructuredMesh2D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "unstructured_2d_dgsem")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "UnstructuredMesh2D" begin
#! format: noindent

@testset "Shallow Water Wet/Dry" begin
    @trixi_testset "elixir_shallowwater_ec.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_ec.jl"),
                            l2=[
                                0.6107326269462766,
                                0.48666631722018877,
                                0.48309775159067053,
                                0.29467422718511704,
                            ],
                            linf=[
                                2.776782342826098,
                                3.2158378644333707,
                                3.652920889487258,
                                2.052861364219655,
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
                                1.2164292510839076,
                                2.6118925543469468e-12,
                                1.1636046671473883e-12,
                                1.2164292510839079,
                            ],
                            linf=[
                                1.5138512282315846,
                                4.998482888288039e-11,
                                2.0246214978154587e-11,
                                1.513851228231574,
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                1.2164292510839085,
                                1.2643106818778908e-12,
                                7.46884905098358e-13,
                                1.2164292510839079,
                            ],
                            linf=[
                                1.513851228231562,
                                1.6287765844373185e-11,
                                6.8766999132716964e-12,
                                1.513851228231574,
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(flux_lax_friedrichs,
                                                                        hydrostatic_reconstruction_audusse_etal),
                                          flux_nonconservative_audusse_etal),
                            tspan=(0.0, 0.2),
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with flux_nonconservative_ersing_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                1.2164292510839083,
                                2.590643638636187e-12,
                                1.0945471514840143e-12,
                                1.2164292510839079,
                            ],
                            linf=[
                                1.5138512282315792,
                                5.0276441977281156e-11,
                                1.9816934589292803e-11,
                                1.513851228231574,
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_ersing_etal),
                            volume_flux=(flux_wintermeyer_etal,
                                         flux_nonconservative_ersing_etal),
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

    @trixi_testset "elixir_shallowwater_source_terms.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.001118134082248467,
                                0.044560486817464634,
                                0.01430926600634214,
                                5.089218476759981e-6,
                            ],
                            linf=[
                                0.007798727223654822,
                                0.34782952734839157,
                                0.11161614702628064,
                                2.6407324614341476e-5,
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.001119719074048042,
                                0.015430106961377458,
                                0.017081080392855975,
                                5.0892184767572974e-6,
                            ],
                            linf=[
                                0.014300786288935274,
                                0.12782194018738569,
                                0.1762517408268396,
                                2.640732461456352e-5,
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_nonconservative_ersing_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.001118046975499805,
                                0.04455969246244461,
                                0.014298120235633432,
                                5.089218476759981e-6,
                            ],
                            linf=[
                                0.007776521213640031,
                                0.34768318303226353,
                                0.11075311228066198,
                                2.6407324614341476e-5,
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_hll" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.0011197190740480426,
                                0.015430106961377272,
                                0.017081080392856225,
                                5.0892184767572974e-6,
                            ],
                            linf=[
                                0.014300786288934386,
                                0.1278219401874079,
                                0.1762517408268396,
                                2.640732461456352e-5,
                            ],
                            surface_flux=(FluxHLL(min_max_speed_naive),
                                          flux_nonconservative_fjordholm_etal),
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

    @trixi_testset "elixir_shallowwater_dirichlet.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_dirichlet.jl"),
                            l2=[
                                1.1577518608940115e-5,
                                4.867189932537344e-13,
                                4.647273240470541e-13,
                                1.1577518608933468e-5,
                            ],
                            linf=[
                                8.394063878602864e-5,
                                1.1469760027632646e-10,
                                1.1146619484429974e-10,
                                8.394063879602065e-5,
                            ],
                            tspan=(0.0, 2.0),
                            surface_flux=(FluxHLL(min_max_speed_naive),
                                          flux_nonconservative_fjordholm_etal),
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

    @trixi_testset "elixir_shallowwater_wall_bc_shockcapturing.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_wall_bc_shockcapturing.jl"),
                            l2=[
                                0.04444388691670699,
                                0.1527771788033111,
                                0.1593763537203512,
                                6.225080476986749e-8,
                            ],
                            linf=[
                                0.6526506870169639,
                                1.980765893182952,
                                2.4807635459119757,
                                3.982097158683473e-7,
                            ],
                            tspan=(0.0, 0.05),
                            surface_flux=(FluxHydrostaticReconstruction(FluxHLL(min_max_speed_naive),
                                                                        hydrostatic_reconstruction_audusse_etal),
                                          flux_nonconservative_audusse_etal))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        let
            t = sol.t[end]
            u_ode = sol.u[end]
            du_ode = similar(u_ode)
            @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
        end
    end

    @trixi_testset "elixir_shallowwater_ec_shockcapturing.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_ec_shockcapturing.jl"),
                            l2=[
                                0.612551520607341,
                                0.5039173660221961,
                                0.49136517934903523,
                                0.29467422718511704,
                            ],
                            linf=[
                                2.7636771472622197,
                                3.236168963021072,
                                3.3363936775653826,
                                2.052861364219655,
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

    @trixi_testset "elixir_shallowwater_three_mound_dam_break.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_three_mound_dam_break.jl"),
                            l2=[
                                0.0892957892027502,
                                0.30648836484407915,
                                2.28712547616214e-15,
                                0.0008778654298684622,
                            ],
                            linf=[
                                0.850329472915091,
                                2.330631694956507,
                                5.783660020252348e-14,
                                0.04326237921249021,
                            ],
                            basis=LobattoLegendreBasis(3),
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
                            l2=[0.0007935561625451243, 0.008825315509943844,
                                0.002429969315645897,
                                0.0007580145888686304, 0.004495741879625235,
                                0.0015758146898767814,
                                6.849532064729749e-6],
                            linf=[0.0059205195991136605, 0.08072126590166251,
                                0.03463806075399023,
                                0.005884818649227186, 0.042658506561995546,
                                0.014125956138838602, 2.5829318284764646e-5],
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
                            l2=[4.706532184998499e-16, 1.1215950712872183e-15,
                                6.7822712922421565e-16,
                                0.002192812926266047, 5.506855295923691e-15,
                                3.3105180099689275e-15,
                                0.0021928129262660085],
                            linf=[4.468647674116255e-15, 1.3607872120431166e-14,
                                9.557155049520056e-15,
                                0.024280130945632084, 6.68910907640583e-14,
                                4.7000983997100496e-14,
                                0.024280130945632732],
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
                            l2=[0.012447632879122346, 0.012361250464676683,
                                0.0009551519536340908,
                                0.09119400061322577, 0.015276216721920347,
                                0.0012126995108983853, 0.09991983966647647],
                            linf=[0.044305765721807444, 0.03279620980615845,
                                0.010754320388190101,
                                0.111309922939555, 0.03663360204931427,
                                0.014332822306649284,
                                0.10000000000000003],
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
                                0.0004897160984647814,
                                0.0005096737100476838,
                                0.000386126716413923,
                                0.0037553984372061602,
                                0.0008470734515158057,
                                0.000914391804695984,
                                0.0039026859922448687,
                                0.0008458945196610379,
                                0.0009797369075336638,
                                1.4701175107332196e-5,
                            ],
                            linf=[
                                0.003156020123921799,
                                0.0023782353180252236,
                                0.002165877463488175,
                                0.027549055479562767,
                                0.007426673513840409,
                                0.0057007815151305374,
                                0.03394812114803547,
                                0.00528271152947668,
                                0.008033142366791368,
                                5.280913235339302e-5,
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
                                0.0002738452813944285,
                                0.00024728982452184266,
                                0.00018186742817211048,
                                0.000515064164935216,
                                0.00020716566255560428,
                                0.0001967259461755342,
                                0.0006765281760969105,
                                0.00026967321034573937,
                                0.0002097997387723772,
                                1.4701175107334422e-5,
                            ],
                            linf=[
                                0.0016143270108415209,
                                0.0013659806938262076,
                                0.0008667976189581372,
                                0.0024878154391543283,
                                0.0007750025239923741,
                                0.0012310851890787178,
                                0.0031150595992053276,
                                0.0015746347140027095,
                                0.0011906648185069368,
                                5.280913235350404e-5,
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced.jl"),
                            l2=[
                                9.267601508692614e-18,
                                1.3079746419456341e-17,
                                1.1469870592306251e-17,
                                0.0021928129262660054,
                                5.551772799657207e-18,
                                5.314204778937538e-18,
                                5.555352319653062e-18,
                                4.268241928252567e-17,
                                6.396791646368386e-18,
                                6.334121301133576e-18,
                                6.427996104426942e-18,
                                5.40279793335249e-17,
                                0.0021928129262660024,
                            ],
                            linf=[
                                2.220446049250313e-16,
                                2.0122792321330962e-16,
                                2.2898349882893854e-16,
                                0.024280130945632805,
                                6.770402623811997e-17,
                                6.668606509370906e-17,
                                6.930883046118016e-17,
                                8.016754559973626e-16,
                                1.2223596839836688e-16,
                                1.0879445133597863e-16,
                                8.534611214117421e-17,
                                9.15698049490698e-16,
                                0.024280130945632836,
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced.jl with flux_lax_friedrichs" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced.jl"),
                            l2=[
                                7.747573306861602e-18,
                                1.217768276787359e-17,
                                9.869557273479402e-18,
                                0.002192812926266005,
                                5.611298116340555e-18,
                                6.0407287687973305e-18,
                                5.968796919323908e-18,
                                4.954396728449881e-17,
                                5.4141534116082716e-18,
                                5.946309178467258e-18,
                                6.157681863345819e-18,
                                4.887115264719831e-17,
                                0.0021928129262660024,
                            ],
                            linf=[
                                3.469446951953614e-17,
                                4.163336342344337e-17,
                                1.3183898417423734e-16,
                                0.02428013094563286,
                                5.852343025496837e-17,
                                6.463560282263022e-17,
                                7.349978245944825e-17,
                                5.309959794423629e-16,
                                4.541940389299994e-17,
                                6.658319696825282e-17,
                                9.223798706573428e-17,
                                5.425491285882212e-16,
                                0.024280130945632836,
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break.jl"),
                            l2=[
                                0.008133096076522545,
                                0.009059379390936938,
                                0.0156472431836693,
                                0.006787446774526842,
                                0.0067495775839911685,
                                0.007737394569507473,
                                0.0001835435993951076,
                                0.00018347939156198144,
                                0.00023500246849234892,
                                0.004003203849568449,
                            ],
                            linf=[
                                0.0268742067416933,
                                0.05243441379266703,
                                0.1590063807340898,
                                0.028581575095277572,
                                0.024486565205936787,
                                0.03953111184090555,
                                0.002625120329243771,
                                0.0026387938388929785,
                                0.0038007822188701103,
                                0.10000000026183736,
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break.jl with flux_lax_friedrichs" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_dam_break.jl"),
                            l2=[
                                0.0089213866827554,
                                0.009502261925352482,
                                0.017690193299797062,
                                0.006845800468099525,
                                0.006849466485572439,
                                0.00761550757687342,
                                0.0004491107702429648,
                                0.00043570105077582146,
                                0.0005832036370132536,
                                0.004003203849568449,
                            ],
                            linf=[
                                0.04276512135712615,
                                0.054061877851283885,
                                0.18896205511083025,
                                0.018902398329971263,
                                0.017901194463979406,
                                0.025135202096835247,
                                0.004065137704454743,
                                0.0038290081477974037,
                                0.005930216968236107,
                                0.10000000026183736,
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
end # MLSWE
end # UnstructuredMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
