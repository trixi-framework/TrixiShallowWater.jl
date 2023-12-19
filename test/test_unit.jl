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

        # test tuple args
        cons_vars = Trixi.prim2cons((H, v, b), equations)
        entropy_vars = Trixi.cons2entropy(cons_vars, equations)
        @test cons_vars ≈ Trixi.entropy2cons(entropy_vars, equations)
    end
end

# TODO: Probably can be removed. Since we dispatch to ShallowWaterEquations1D we don't need 
#       equation specific wave speed estimates.
# @timed_testset "Connectivity with Trixi.jl" begin
#     u = SVector(SVector(1, 0.5, 0.0))
#     orientation = 1
#     equations = ShallowWaterEquationsWetDry1D(gravity_constant = 9.81)
#     equations_trixi = ShallowWaterEquations1D(gravity_constant = 9.81)

#     # We only need to check equivalence between the equation systems. The functionality is tested in
#     # Trixi.jl. We choose these specific ones to improve code coverage, as they are not tested in
#     # any elixirs.
#     @test min_max_speed_einfeldt(u, u, orientation, equations) ==
#           min_max_speed_einfeldt(u, u, orientation, equations_trixi)
#     @test min_max_speed_davis(u, u, orientation, equations) ==
#           min_max_speed_davis(u, u, orientation, equations_trixi)
#     @test max_abs_speed_naive(u, u, orientation, equations) ==
#           max_abs_speed_naive(u, u, orientation, equations_trixi)
# end
end

end # module
