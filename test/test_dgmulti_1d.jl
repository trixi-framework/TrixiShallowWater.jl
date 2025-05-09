module TestExamplesDGMulti1D

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

EXAMPLES_DIR = pkgdir(TrixiShallowWater, "examples", "dgmulti_1d")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "DGMulti 1D" begin
#! format: noindent

@trixi_testset "elixir_shallow_water_quasi_1d.jl (SBP) " begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "elixir_shallow_water_quasi_1d.jl"),
                        cells_per_dimension=(8,),
                        approximation_type=SBP(),
                        l2=[
                            3.0300196635805022e-6,
                            1.6921833812545857e-5,
                            2.9844594164368975e-16,
                            1.1012004949980629e-15
                        ],
                        linf=[
                            1.2043309307818717e-5,
                            5.346754311919e-5,
                            9.43689570931383e-16,
                            2.220446049250313e-15
                        ])
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    let
        t = sol.t[end]
        u_ode = sol.u[end]
        du_ode = similar(u_ode)
        @test (@allocated Trixi.rhs!(du_ode, u_ode, semi, t)) < 1000
    end
end
end # DGMulti1D

end # module
