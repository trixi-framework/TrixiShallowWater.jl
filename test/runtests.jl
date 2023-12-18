using TrixiShallowWater
using Test

# We run tests in parallel with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
# By default, we just run the threaded tests since they are relatively cheap
# and test a good amount of different functionality.
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")
const TRIXI_MPI_NPROCS = clamp(Sys.CPU_THREADS, 2, 3)
const TRIXI_NTHREADS = clamp(Sys.CPU_THREADS, 2, 3)

@time @testset "TrixiShallowWater.jl tests" begin
    @time if TRIXI_TEST == "all" || TRIXI_TEST == "tree_part1"
        include("test_tree_1d_shallowwater_wet_dry.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "upstream"
        @testset "baz()" begin
            @test TrixiShallowWater.baz() isa String
        end
    end
end
