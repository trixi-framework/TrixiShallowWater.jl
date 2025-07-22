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
                                0.29467422718511704
                            ],
                            linf=[
                                2.776782342826098,
                                3.2158378644333707,
                                3.652920889487258,
                                2.052861364219655
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

    @trixi_testset "elixir_shallowwater_ec_float32.jl" begin
        # Expected errors are nearly all taken from elixir_shallowwater_ec.jl
        @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallowwater_ec_float32.jl"),
                            l2=[
                                Float32(0.6107326269462766),
                                Float32(0.48666631722018877),
                                Float32(0.48309775159067053),
                                Float32(0.29467422718511704)
                            ],
                            linf=[
                                Float32(2.776782342826098),
                                3.2116454f0, # this needs to be adapted
                                3.6616623f0, # this needed to be adapted
                                Float32(2.052861364219655)
                            ],
                            tspan=(0.0f0, 0.25f0),
                            RealT=Float32)
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
                                1.2164292510839079
                            ],
                            linf=[
                                1.5138512282315846,
                                4.998482888288039e-11,
                                2.0246214978154587e-11,
                                1.513851228231574
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
                                1.2164292510839079
                            ],
                            linf=[
                                1.513851228231562,
                                1.6287765844373185e-11,
                                6.8766999132716964e-12,
                                1.513851228231574
                            ],
                            surface_flux=(FluxHydrostaticReconstruction(FluxLaxFriedrichs(max_abs_speed_naive),
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

    @trixi_testset "elixir_shallowwater_well_balanced.jl with flux_nonconservative_wintermeyer_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_well_balanced.jl"),
                            l2=[
                                1.2164292510839083,
                                2.590643638636187e-12,
                                1.0945471514840143e-12,
                                1.2164292510839079
                            ],
                            linf=[
                                1.5138512282315792,
                                5.0276441977281156e-11,
                                1.9816934589292803e-11,
                                1.513851228231574
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_wintermeyer_etal),
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
                                5.089218476759981e-6
                            ],
                            linf=[
                                0.007798727223654822,
                                0.34782952734839157,
                                0.11161614702628064,
                                2.6407324614341476e-5
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
                                5.0892184767572974e-6
                            ],
                            linf=[
                                0.014300786288935274,
                                0.12782194018738569,
                                0.1762517408268396,
                                2.640732461456352e-5
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

    @trixi_testset "elixir_shallowwater_source_terms.jl with flux_nonconservative_wintermeyer_etal" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_source_terms.jl"),
                            l2=[
                                0.001118046975499805,
                                0.04455969246244461,
                                0.014298120235633432,
                                5.089218476759981e-6
                            ],
                            linf=[
                                0.007776521213640031,
                                0.34768318303226353,
                                0.11075311228066198,
                                2.6407324614341476e-5
                            ],
                            surface_flux=(flux_wintermeyer_etal,
                                          flux_nonconservative_wintermeyer_etal),
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
                                5.0892184767572974e-6
                            ],
                            linf=[
                                0.014300786288934386,
                                0.1278219401874079,
                                0.1762517408268396,
                                2.640732461456352e-5
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
                                1.1577518608928104e-5,
                                4.761468345949485e-13,
                                4.546054642605431e-13,
                                1.157751860893347e-5
                            ],
                            linf=[
                                8.394063878847113e-5,
                                1.1211499950000422e-10,
                                1.0890394254534975e-10,
                                8.394063879602065e-5
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
                                6.225080476986749e-8
                            ],
                            linf=[
                                0.6526506870169639,
                                1.980765893182952,
                                2.4807635459119757,
                                3.982097158683473e-7
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
                                0.29467422718511704
                            ],
                            linf=[
                                2.7636771472622197,
                                3.236168963021072,
                                3.3363936775653826,
                                2.052861364219655
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
                                0.08930075417161291,
                                0.3065009462584156,
                                1.9814680139136418e-16,
                                0.0008778654298685837
                            ],
                            linf=[
                                0.8504438370731531,
                                2.3307422133452445,
                                3.762330934103131e-15,
                                0.04326237921249021
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

    @trixi_testset "elixir_shallowwater_inflow_outflow.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_inflow_outflow.jl"),
                            l2=[
                                0.05395240720695749,
                                0.17993463107947336,
                                0.17890511959110378,
                                0.0
                            ],
                            linf=[
                                0.10340479800696867,
                                0.3998099631746912,
                                0.4041504686758526,
                                0.0
                            ],
                            tspan=(0.0, 3.0))
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

    @trixi_testset "elixir_shallowwater_twolayer_dam_break.jl with FluxLaxFriedrichs(max_abs_speed_naive)" begin
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
                                2.845107067734722e-5,
                                2.541332828758267e-5,
                                1.3776844294886904e-5,
                                3.21641609844227e-5,
                                1.9645259635774905e-5,
                                1.2494077534998188e-5,
                                3.842168759551814e-5,
                                2.3940963660476904e-5,
                                1.5539403052020932e-5,
                                7.244891628173627e-7
                            ],
                            linf=[
                                0.0001686597682639679,
                                0.00016943475945224717,
                                0.00010161333860780886,
                                0.00018234203945666216,
                                0.0001440335843042595,
                                7.821014430642315e-5,
                                0.0002159224957090089,
                                0.00016273539529820802,
                                9.798918563769243e-5,
                                4.040896422807805e-6
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
                                1.918097575328871e-5,
                                1.502019981850744e-5,
                                5.042301180852648e-6,
                                1.9602948809781437e-5,
                                1.1406394406547229e-5,
                                5.066224008439183e-6,
                                2.3354649343891286e-5,
                                1.433469494607502e-5,
                                6.046742009426357e-6,
                                7.244891628173627e-7
                            ],
                            linf=[
                                0.0001267907408812885,
                                0.00013376989668978378,
                                3.3611013632417475e-5,
                                0.00013760680376173617,
                                0.00011659018473553218,
                                3.15196814599239e-5,
                                0.00017067740318088553,
                                0.00012830136581865048,
                                3.6074428482746335e-5,
                                4.040896422807805e-6
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

    @trixi_testset "elixir_shallowwater_multilayer_convergence.jl with FluxHydrostaticReconstruction" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence.jl"),
                            l2=[
                                2.8451070677919015e-5,
                                2.5413328287320148e-5,
                                1.3776844295056111e-5,
                                3.216416098536785e-5,
                                1.964525963539178e-5,
                                1.2494077535320725e-5,
                                3.842168759685933e-5,
                                2.3940963660143302e-5,
                                1.553940305228399e-5,
                                7.244891628173627e-7
                            ],
                            linf=[
                                0.00016865976829461005,
                                0.0001694347594369816,
                                0.00010161333860081445,
                                0.00018234203943068295,
                                0.00014403358429010416,
                                7.821014428416317e-5,
                                0.0002159224956841399,
                                0.00016273539534478187,
                                9.798918562298198e-5,
                                4.040896422807805e-6
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
                                0.0021928129262660024
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
                                0.024280130945632836
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
                                0.0021928129262660024
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
                                0.024280130945632836
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

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced_wet_dry.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced_wet_dry.jl"),
                            l2=[
                                0.258337877456028,
                                0.4536334309305916,
                                0.2424975342044823,
                                1.861822857108904e-15,
                                1.0728615421620665e-15,
                                4.893078129679801e-16,
                                2.951473837608095e-15,
                                9.819198357343179e-16,
                                3.666055853677385e-16,
                                0.9049896638396252
                            ],
                            linf=[
                                0.5000000000000003,
                                0.5000000000000064,
                                0.485509867641643,
                                1.7168360326077987e-14,
                                6.4508928741529934e-15,
                                3.663826180637297e-15,
                                2.370261384584603e-14,
                                6.3907606640361836e-15,
                                2.0890952713058902e-15,
                                1.3458935664973586
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
                                0.008133096076522545,
                                0.009059379390936938,
                                0.0156472431836693,
                                0.006787446774526842,
                                0.0067495775839911685,
                                0.007737394569507473,
                                0.0001835435993951076,
                                0.00018347939156198144,
                                0.00023500246849234892,
                                0.004003203849568449
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
                                0.10000000026183736
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

    @trixi_testset "elixir_shallowwater_multilayer_dam_break.jl with FluxLaxFriedrichs(max_abs_speed_naive)" begin
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
                                0.004003203849568449
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
                                0.10000000026183736
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
                                0.025056397690317766,
                                0.02413438932652429,
                                0.05305193031437963,
                                0.01570655272565154,
                                0.015434075201987807,
                                0.03802395905825049,
                                0.0035119966268674424,
                                0.003140567605820363,
                                0.0047903431628184764,
                                0.00400367487891041
                            ],
                            linf=[
                                0.13756514254828095,
                                0.13141572813228075,
                                0.40279814892528143,
                                0.05125099836463482,
                                0.05082705371717001,
                                0.20100787456634925,
                                0.021397757474313436,
                                0.019301532073029846,
                                0.031247259773652267,
                                0.10003205938749304
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
end # UnstructuredMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
