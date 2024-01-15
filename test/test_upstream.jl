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
    @trixi_testset "1D-Test: elixir_shallowwater_well_balanced_nonperiodic.jl with wall boundary" begin
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
end

# Clean up afterwards: delete output directory
@test_nowarn rm(outdir, recursive = true)
end # Upstream tests

end # module
