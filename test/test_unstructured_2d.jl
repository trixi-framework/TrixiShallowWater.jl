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
                                     l2 = [1.8656084144845292e-5, 1.917796508073001e-5, 1.7454512418793462e-5, 0.00014048665094953486, 3.1328989570708975e-5, 3.512474430566952e-5, 0.00014804378967592363, 3.2299027874400884e-5, 3.881063343808959e-5, 2.393199280374949e-7],
                                     linf = [0.00013984354089258133, 0.00011014730412189921, 0.00010077985342477058, 0.0011054081397450233, 0.00022137639447733504, 0.0002451045776211691, 0.0009651462415440903, 0.00026431961702122475, 0.0002728818589384785, 1.03643740911874e-6],
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
                                     l2 = [1.0098144596783888e-5, 8.667598940295975e-6, 7.5615805463112364e-6, 1.8150527998748446e-5, 6.935960840543087e-6, 7.335149149911856e-6, 2.403682140927019e-5, 9.684693895783513e-6, 8.27535107484968e-6, 2.393199280374949e-7],
                                     linf = [6.62853303028399e-5, 5.144519565108974e-5, 3.6540170641030656e-5, 0.00013228821066646468, 4.67241008220709e-5, 4.623162831440819e-5, 0.00015925845769571012, 7.204370895297352e-5, 4.566973688335807e-5, 1.03643740911874e-6],
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
                                     l2 = [0.008264086240458098, 0.008904024684630463, 0.015040884852261934, 0.0068380853981447965, 0.006803190764735809, 0.007681935438533932, 0.00014126104190197007, 0.00013800738645477302, 0.00015518676788954237, 0.0062500008375575705],
                                     linf = [0.0271281524526652, 0.038446397522241216, 0.10189873392468171, 0.02279784901518067, 0.02344053175169211, 0.0384823318069141, 0.0015613916350439184, 0.0015293762865072726, 0.0020984165518842324, 0.05808926980412543],
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
                                     l2 = [0.008289485483638842, 0.008744140233824497, 0.014162544835509572, 0.0067108677626928565, 0.006699962946263383, 0.0075441143283640315, 0.00012772255623776208, 0.0001243555119128215, 0.00013480820924093037, 0.0062500008375575705],
                                     linf = [0.024574998802901732, 0.024386474982493633, 0.09088941363545519, 0.01828534534858913, 0.017626892238781267, 0.02546243728251761, 0.0007231987972067275, 0.0007131387402676157, 0.0007959906366350036, 0.05808926980412543],
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
