using Trixi
using TrixiShallowWater
using Test

# We run tests with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")

@time @testset "TrixiShallowWater.jl tests" begin
    @time if TRIXI_TEST == "all"
        include("test_tree_1d.jl")
        include("test_tree_2d.jl")
        include("test_unstructured_2d.jl")
        include("test_structured_2d.jl")
        include("test_unit.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "upstream"
        @testset "Namespace conflicts" begin
            # Test for namespace conflicts between TrixiShallowWater.jl and Trixi.jl
            for name in names(Trixi)
                @test !(name in names(TrixiShallowWater))
            end
        end

        # Run upstream tests for each mesh and dimension to test compatibility with Trixi.jl
        include("test_upstream.jl")
    end
end
