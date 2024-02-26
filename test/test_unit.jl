module TestUnit

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

# Start with a clean environment: remove TrixiShallowWater.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

# Run various unit (= non-elixir-triggered) tests
@testset "Unit tests" begin
#! format: noindent

@timed_testset "Printing indicators/controllers" begin
    # OBS! Constructing indicators/controllers using the parameters below doesn't make sense. It's
    # just useful to run basic tests of `show` methods.

    indicator_hg_swe = IndicatorHennemannGassnerShallowWater(1.0, 0.0, true, "variable",
                                                             "cache")
    @test_nowarn show(stdout, indicator_hg_swe)
end

@timed_testset "Shallow water conversion between conservative/entropy variables" begin
    H, v1, v2, b = 3.5, 0.25, 0.1, 0.4

    let equations = ShallowWaterEquationsWetDry1D(gravity_constant = 9.8)
        cons_vars = prim2cons(SVector(H, v1, b), equations)
        entropy_vars = cons2entropy(cons_vars, equations)
        @test cons_vars ≈ entropy2cons(entropy_vars, equations)

        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
        # test tuple args
        cons_vars = prim2cons((H, v1, b), equations)
        entropy_vars = cons2entropy(cons_vars, equations)
        @test cons_vars ≈ entropy2cons(entropy_vars, equations)
    end

    let equations = ShallowWaterEquationsWetDry2D(gravity_constant = 9.8)
        cons_vars = prim2cons(SVector(H, v1, v2, b), equations)
        entropy_vars = cons2entropy(cons_vars, equations)
        @test cons_vars ≈ entropy2cons(entropy_vars, equations)

        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)

        # test tuple args
        cons_vars = prim2cons((H, v1, v2, b), equations)
        entropy_vars = cons2entropy(cons_vars, equations)
        @test cons_vars ≈ entropy2cons(entropy_vars, equations)
    end
end

@timed_testset "Consistency check for waterheight_pressure" begin
    H, v1, v2, b = 3.5, 0.25, 0.1, 0.4

    let equations = ShallowWaterEquationsWetDry1D(gravity_constant = 9.8)
        cons_vars = prim2cons(SVector(H, v1, b), equations)
        @test waterheight_pressure(cons_vars, equations) ≈
              Trixi.waterheight(cons_vars, equations) * pressure(cons_vars, equations)
    end

    let equations = ShallowWaterEquationsWetDry2D(gravity_constant = 9.8)
        cons_vars = prim2cons(SVector(H, v1, v2, b), equations)
        @test waterheight_pressure(cons_vars, equations) ≈
              Trixi.waterheight(cons_vars, equations) * pressure(cons_vars, equations)
    end
end
end

end # module
