module TestUnit

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

# Run various unit (= non-elixir-triggered) tests
@testset "Unit tests" begin
#! format: noindent

@timed_testset "Shallow water conversion between conservative/entropy variables" begin
    H, v, b = 3.5, 0.25, 0.4

    let equations = ShallowWaterEquationsWetDry1D(gravity_constant = 9.8)
        cons_vars = Trixi.prim2cons(SVector(H, v, b), equations)
        entropy_vars = Trixi.cons2entropy(cons_vars, equations)
        @test cons_vars ≈ Trixi.entropy2cons(entropy_vars, equations)

        total_energy = Trixi.energy_total(cons_vars, equations)
        @test total_energy ≈ Trixi.entropy(cons_vars, equations)
        @test total_energy ≈
              Trixi.energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
        # test tuple args
        cons_vars = Trixi.prim2cons((H, v, b), equations)
        entropy_vars = Trixi.cons2entropy(cons_vars, equations)
        @test cons_vars ≈ Trixi.entropy2cons(entropy_vars, equations)
    end
end
end

end # module
