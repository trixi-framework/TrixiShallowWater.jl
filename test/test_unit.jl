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

@timed_testset "Two-layer shallow water conversion between conservative/entropy variables" begin
    H_upper, v1_upper, v2_upper, H_lower, v1_lower, v2_lower, b = 3.5, 0.25, 0.13, 2.5,
                                                                  0.1, 0.37, 0.4

    let equations = ShallowWaterTwoLayerEquations1D(gravity_constant = 9.8,
                                                    rho_upper = 0.9, rho_lower = 1.0)
        cons_vars = prim2cons(SVector(H_upper, v1_upper, H_lower, v1_lower, b),
                              equations)
        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
    end

    let equations = ShallowWaterTwoLayerEquations2D(gravity_constant = 9.8,
                                                    rho_upper = 0.9, rho_lower = 1.0)
        cons_vars = prim2cons(SVector(H_upper, v1_upper, v2_upper, H_lower, v1_lower,
                                      v2_lower, b), equations)
        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
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

# Test consistency of the wrapper functions for the wave speed estimates
@timed_testset "Consistency check for wave speed estimates of the SWE" begin
    let equations = ShallowWaterEquationsWetDry1D(gravity_constant = 9.8)
        u_rr = SVector(1.2, 0.3, 0.7)
        u_ll = SVector(0.25, 0.1, 0.4)
        orientation = 1

        @test min_max_speed_naive(u_ll, u_rr, orientation,
                                  equations) ==
              min_max_speed_naive(u_ll, u_rr, orientation,
                                  equations.basic_swe)
        @test min_max_speed_davis(u_ll, u_rr, orientation,
                                  equations) ==
              min_max_speed_davis(u_ll, u_rr, orientation,
                                  equations.basic_swe)
        @test min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                     equations) ==
              min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                     equations.basic_swe)
        @test min_max_speed_naive(u_ll, u_rr, orientation,
                                  equations) ==
              min_max_speed_naive(u_ll, u_rr, orientation,
                                  equations.basic_swe)
    end

    let equations = ShallowWaterEquationsWetDry2D(gravity_constant = 9.8)
        u_rr = SVector(1.2, 0.3, 0.2, 0.7)
        u_ll = SVector(0.25, 0.1, 0.3, 0.4)
        orientations = [1, 2]
        normal_directions = [SVector(1.0, 0.0),
            SVector(0.0, 1.0),
            SVector(0.5, -0.5),
            SVector(-1.2, 0.3)]

        for orientation in orientations
            @test min_max_speed_naive(u_ll, u_rr, orientation,
                                      equations) ==
                  min_max_speed_naive(u_ll, u_rr, orientation,
                                      equations.basic_swe)
            @test min_max_speed_davis(u_ll, u_rr, orientation,
                                      equations) ==
                  min_max_speed_davis(u_ll, u_rr, orientation,
                                      equations.basic_swe)
            @test min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                         equations) ==
                  min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                         equations.basic_swe)
            @test min_max_speed_naive(u_ll, u_rr, orientation,
                                      equations) ==
                  min_max_speed_naive(u_ll, u_rr, orientation,
                                      equations.basic_swe)
        end

        for normal_direction in normal_directions
            @test min_max_speed_naive(u_ll, u_rr, normal_direction,
                                      equations) ==
                  min_max_speed_naive(u_ll, u_rr, normal_direction,
                                      equations.basic_swe)
            @test min_max_speed_davis(u_ll, u_rr, normal_direction,
                                      equations) ==
                  min_max_speed_davis(u_ll, u_rr, normal_direction,
                                      equations.basic_swe)
            @test min_max_speed_einfeldt(u_ll, u_rr, normal_direction,
                                         equations) ==
                  min_max_speed_einfeldt(u_ll, u_rr, normal_direction,
                                         equations.basic_swe)
            @test min_max_speed_naive(u_ll, u_rr, normal_direction,
                                      equations) ==
                  min_max_speed_naive(u_ll, u_rr, normal_direction,
                                      equations.basic_swe)
        end
    end
end

@timed_testset "Exception check for the 2LSWE" begin
    error_message = "Invalid input: Densities must be chosen such that rho_upper < rho_lower"
    @test_throws error_message ShallowWaterTwoLayerEquations1D(gravity_constant = 9.81,
                                                               rho_upper = 1.0,
                                                               rho_lower = 0.9)
    @test_throws error_message ShallowWaterTwoLayerEquations2D(gravity_constant = 9.81,
                                                               rho_upper = 1.0,
                                                               rho_lower = 0.9)
end

@timed_testset "Input argument check for the ML-SWE" begin
    @test_throws ArgumentError ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                                                 rhos = [
                                                                     -1.0,
                                                                     0.1,
                                                                     0.2,
                                                                 ])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                                                 rhos = [0.1, 0.3, 0.2])
    # Ensure that both tuple and array input are equivalent    
    @test ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                            rhos = [0.1, 0.2, 0.3]) ==
          ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                            rhos = (0.1, 0.2, 0.3))
end
end # Unit tests

end # module
