using Trixi
using TrixiShallowWater
using Test
using MPI: mpiexec

# We run tests with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")
const TRIXI_MPI_NPROCS = clamp(Sys.CPU_THREADS, 2, 3)

@time @testset "TrixiShallowWater.jl tests" begin
    # This is placed first since tests error out otherwise if `TRIXI_TEST == "all"`,
    # at least on some systems.
    @time if TRIXI_TEST == "all" || TRIXI_TEST == "mpi"
        # Do a dummy `@test true`:
        # If the process errors out the testset would error out as well,
        # cf. https://github.com/JuliaParallel/MPI.jl/pull/391
        @test true

        # We provide a `--heap-size-hint` to avoid/reduce out-of-memory errors during CI testing
        mpiexec() do cmd
            run(`$cmd -n $TRIXI_MPI_NPROCS $(Base.julia_cmd()) --threads=1 --check-bounds=yes --heap-size-hint=0.5G $(abspath("test_mpi.jl"))`)
        end
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "tree_1d"
        include("test_tree_1d.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "tree_2d"
        include("test_tree_2d.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "unstructured_2d"
        include("test_unstructured_2d.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "structured_2d"
        include("test_structured_2d.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "p4est_2d"
        include("test_p4est_2d.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "unit"
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
