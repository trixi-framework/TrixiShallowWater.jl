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
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
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
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
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
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
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
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    # Note, these values may change as the functionality of well-balanced mortars
    # with AMR and wet/dry are further developed according to the issue
    # https://github.com/trixi-framework/TrixiShallowWater.jl/issues/77
    # This is run longer to activate the pieces of the refine portion of the limiter
    @trixi_testset "elixir_shallowwater_perturbation_wet_dry_amr.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_perturbation_wet_dry_amr.jl"),
                            l2=[
                                0.38992332005657065,
                                0.2658037399348242,
                                0.27166206980738195,
                                0.45643510424305106
                            ],
                            linf=[
                                2.030377008431841,
                                1.930801096724945,
                                2.173823340700309,
                                0.7495177590247986
                            ],
                            tspan=(0.0, 0.121))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    @trixi_testset "elixir_shallowwater_rainfall_inclined_plane.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_rainfall_inclined_plane.jl"),
                            l2=[
                                0.00032207028548645125,
                                2.30277557795753e-6,
                                4.44997187575836e-19,
                                1.8433475390376417e-15
                            ],
                            linf=[
                                0.000324115924751592,
                                2.3248063411922284e-6,
                                1.9913550038729716e-17,
                                8.881784197001252e-15
                            ],
                            precipitation_rate=(x, t) -> 1e-3 * t,
                            tspan=(0.0, 1.0))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end
end # SWE

@testset "Multilayer Shallow Water" begin
    @trixi_testset "elixir_shallowwater_multilayer_convergence_sc_subcell.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence_sc_subcell.jl"),
                            l2=[
                                0.001004593980644189,
                                0.0006189233371193447,
                                0.000140719252286664,
                                0.0009194923081820151,
                                0.0004928204355037366,
                                0.00013565779779670765,
                                0.0010852948428510846,
                                0.0006186124110474113,
                                0.00015035777111446134,
                                0.00019675440964321886
                            ],
                            linf=[
                                0.005808422635784183,
                                0.00356346714360839,
                                0.0008266733071794485,
                                0.004054746709408086,
                                0.0028567470129725603,
                                0.0005170835013522668,
                                0.004953411660586049,
                                0.0035121047612464706,
                                0.0006376457357819554,
                                0.00043748911723784367
                            ],
                            tspan=(0.0, 0.1))
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        # Larger values for allowed allocations due to usage of custom
        # integrator which are not *recorded* for the methods from
        # OrdinaryDiffEq.jl
        # Corresponding issue: https://github.com/trixi-framework/Trixi.jl/issues/1877
        @test_allocations(Trixi.rhs!, semi, sol, 15000)
    end

    @trixi_testset "elixir_shallowwater_multilayer_convergence_sc_subcell_curved.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_convergence_sc_subcell_curved.jl"),
                            l2=[
                                0.00013275624153302434,
                                0.0001586600356395913,
                                0.000158660035639501,
                                2.9583315272612922e-5
                            ],
                            linf=[
                                0.0007544272991792944,
                                0.0007250877164874936,
                                0.00072508771648927,
                                9.31948456051046e-5
                            ])
        # Allocation testing is disabled as the symbolic source term computation is known to cause
        # allocations.
    end

    @trixi_testset "elixir_shallowwater_multilayer_well_balanced_wet_dry_nonconforming.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_well_balanced_wet_dry_nonconforming.jl"),
                            l2=[
                                0.17389058166483545,
                                1.3442539891317369e-15,
                                1.390606878462557e-15,
                                0.18415049792609314],
                            linf=[
                                0.4160052818017864,
                                1.0894983713455818e-14,
                                1.132829841145215e-14,
                                0.41600528180178625],
                            tspan=(0.0, 0.25),
                            atol=1e-10)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    # Note, these values may change as the functionality of well-balanced mortars
    # with AMR and wet/dry are further developed according to the issue
    # https://github.com/trixi-framework/TrixiShallowWater.jl/issues/77
    @trixi_testset "elixir_shallowwater_multilayer_perturbation_wet_dry_amr.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_perturbation_wet_dry_amr.jl"),
                            l2=[
                                0.4177955194235321,
                                0.2818773677696145,
                                0.33178075151092923,
                                0.45643505879802193],
                            linf=[
                                2.129878820442652,
                                3.292431218004944,
                                3.3701108685790135,
                                0.7495177590247986],
                            tspan=(0.0, 0.01))

        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    # Note, complex example that uses TrixiBottomTopography.jl to approximate
    # bathymetry and wave maker boundary condition. Tests several components
    # of the TrixiShallowWater.jl toolchain. Note, does not run long enough
    # for the AMR to fire.
    @trixi_testset "elixir_shallowwater_multilayer_monai_flood_amr.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_monai_flood_amr.jl"),
                            l2=[
                                0.00030685696894658095,
                                3.869268491914826e-6,
                                6.345699712019463e-11,
                                0.00030875165148307556],
                            linf=[
                                0.0043875707537211345,
                                3.562731345127931e-5,
                                1.1903363999889186e-9,
                                0.004387569562413221],
                            tspan=(0.0, 0.25),
                            atol=5e-12)

        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    @trixi_testset "elixir_shallowwater_multilayer_three_mound_dam_break_amr.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_three_mound_dam_break_amr.jl"),
                            l2=[
                                0.13070149058594405,
                                0.3557070936603237,
                                3.391658168396266e-14,
                                0.001960648080166412],
                            linf=[
                                1.1241019512063912,
                                2.4986126972245994,
                                3.782627199174957e-12,
                                0.044079775502972124],
                            tspan=(0.0, 0.3),
                            atol=1e-10)

        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    @trixi_testset "elixir_shallowwater_multilayer_blast_wet_dry_sc_subcell.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_blast_wet_dry_sc_subcell.jl"),
                            l2=[
                                0.31374291044904945,
                                0.9623170373390161,
                                0.9623163985536942,
                                1.167749829555993e-16
                            ],
                            linf=[
                                1.4347685747786731,
                                4.760470726550078,
                                4.761142120243878,
                                4.440892098500626e-16
                            ],
                            tspan=(0.0, 0.05),
                            # Increase the absolute tolerance to account for varying results with
                            # with the two-sided limiter on different architectures.
                            # See https://github.com/trixi-framework/Trixi.jl/pull/2007
                            atol=5e-4)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        # Larger values for allowed allocations due to usage of custom
        # integrator which are not *recorded* for the methods from
        # OrdinaryDiffEq.jl
        # Corresponding issue: https://github.com/trixi-framework/Trixi.jl/issues/1877
        @test_allocations(Trixi.rhs!, semi, sol, 15000)
    end

    @trixi_testset "elixir_shallowwater_multilayer_blast_wet_dry_sc_subcell with periodic BC.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_blast_wet_dry_sc_subcell.jl"),
                            l2=[
                                0.31374291044904945,
                                0.9623170373390161,
                                0.9623163985536942,
                                1.167749829555993e-16
                            ],
                            linf=[
                                1.4347685747786731,
                                4.760470726550078,
                                4.761142120243878,
                                4.440892098500626e-16
                            ],
                            tspan=(0.0, 0.05),
                            mesh=P4estMesh(trees_per_dimension, polydeg = 3,
                                           coordinates_min = coordinates_min,
                                           coordinates_max = coordinates_max,
                                           initial_refinement_level = 5,
                                           periodicity = true),
                            semi=SemidiscretizationHyperbolic(mesh, equations,
                                                              initial_condition, solver,
                                                              boundary_conditions = boundary_condition_periodic),
                            # Increase the absolute tolerance to account for varying results with
                            # with the two-sided limiter on different architectures.
                            # See https://github.com/trixi-framework/Trixi.jl/pull/2007
                            atol=5e-4)
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        # Larger values for allowed allocations due to usage of custom
        # integrator which are not *recorded* for the methods from
        # OrdinaryDiffEq.jl
        # Corresponding issue: https://github.com/trixi-framework/Trixi.jl/issues/1877
        @test_allocations(Trixi.rhs!, semi, sol, 15000)
    end

    @trixi_testset "elixir_shallowwater_multilayer_ec_sc_subcell.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_multilayer_ec_sc_subcell.jl"),
                            l2=[
                                0.04338581073333642,
                                0.12068863355590975,
                                0.12068863355590971,
                                7.906244739074657e-18
                            ],
                            linf=[
                                0.2550462084344076,
                                0.42004261172048024,
                                0.4200426117204863,
                                1.1102230246251565e-16
                            ])
        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        # Larger values for allowed allocations due to usage of custom
        # integrator which are not *recorded* for the methods from
        # OrdinaryDiffEq.jl
        # Corresponding issue: https://github.com/trixi-framework/Trixi.jl/issues/1877
        @test_allocations(Trixi.rhs!, semi, sol, 15000)
    end
end # MLSWE

@testset "Shallow Water Exner" begin
    @trixi_testset "elixir_shallowwater_exner_convergence_grass.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_exner_convergence_grass.jl"),
                            l2=[
                                0.0007843486672958833,
                                0.006981545959119388,
                                0.007421898782394062,
                                0.00029960154134387536
                            ],
                            linf=[
                                0.007656547860054097,
                                0.0787603790906501,
                                0.10597723748192145,
                                0.0013988745083841625
                            ],
                            tspan=(0.0, 0.1))

        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    @trixi_testset "elixir_shallowwater_exner_ec_mpm.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_exner_ec_mpm.jl"),
                            l2=[
                                0.9861501067398692,
                                0.5932601231318787,
                                0.6914651918642902,
                                0.45633443286763237
                            ],
                            linf=[
                                2.3936795222158076,
                                5.8907408754608515,
                                4.337103147537311,
                                2.9436448052630597
                            ],
                            tspan=(0.0, 0.05))

        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end

    @trixi_testset "elixir_shallowwater_exner_channel_mpm.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_exner_channel_mpm.jl"),
                            l2=[
                                0.021627221304021577,
                                0.0008251875449399632,
                                0.07133184657604669,
                                2.9459714466880464e-7
                            ],
                            linf=[
                                0.1715300794040595,
                                0.05307818092229056,
                                0.57807103428298,
                                1.2915591585352382e-5
                            ],
                            tspan=(0.0, 0.1))

        # Ensure that we do not have excessive memory allocations
        # (e.g., from type instabilities)
        @test_allocations(Trixi.rhs!, semi, sol, 1000)
    end
end # Shallow Water Exner
end # P4estMesh2D

# Clean up afterwards: delete TrixiShallowWater.jl output directory
@test_nowarn rm(outdir, recursive = true)

end # module
