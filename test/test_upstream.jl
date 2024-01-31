module TestExamplesUpstream

using Test
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples")

# Start with a clean environment: remove output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

# Run upstream tests for each mesh and dimension to test compatibility with Trixi.jl
@testset "Upstream tests" begin
#! format: noindent

# Run tests for TreeMesh
@testset "TreeMesh" begin
    # Shallow water wet/dry 1D
    @trixi_testset "TreeMesh1D: elixir_shallowwater_well_balanced_nonperiodic.jl with wall boundary" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR, "tree_1d_dgsem",
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
    # TODO: add upstream tests in 2D and positivity-preserving tests

    # Shallow water wet/dry 2D
    # TreeMesh2D
    @trixi_testset "TreeMesh2D: elixir_shallowwater_conical_island.jl" begin
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
    # Unstructured2D
    @trixi_testset "Unstructured2D: elixir_shallowwater_wall_bc_shockcapturing.jl" begin
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
    # Structured2D
    @trixi_testset "Structured2D: elixir_shallowwater_conical_island.jl" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR,
                                     "elixir_shallowwater_conical_island.jl"),
                            l2=[
                                0.04593154164306353,
                                0.1644534881916908,
                                0.16445348819169076,
                                0.0011537702354532122,
                            ],
                            linf=[
                                0.21100717610846442,
                                0.9501592344310412,
                                0.950159234431041,
                                0.021790250683516296,
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
end

# Clean up afterwards: delete output directory
@test_nowarn rm(outdir, recursive = true)
end # Upstream tests

end # module
