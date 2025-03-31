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
        # Test conversion between primitive and conservative variables
        prim_vars = SVector(H, v1, b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

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
        # Test conversion between primitive and conservative variables
        prim_vars = SVector(H, v1, v2, b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

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
        # Test conversion between primitive and conservative variables
        prim_vars = SVector(H_upper, v1_upper, H_lower, v1_lower, b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
    end

    let equations = ShallowWaterTwoLayerEquations2D(gravity_constant = 9.8,
                                                    rho_upper = 0.9, rho_lower = 1.0)
        # Test conversion between primitive and conservative variables
        prim_vars = SVector(H_upper, v1_upper, v2_upper, H_lower, v1_lower,
                            v2_lower, b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
    end
end

@timed_testset "Multilayer shallow water conversion between conservative/entropy variables" begin
    H = (3.5, 2.5, 1.5)
    v1 = (0.25, 0.1, 0.37)
    v2 = (0.13, 0.2, 0.3)
    b = 0.4

    let equations = ShallowWaterMultiLayerEquations1D(gravity_constant = 9.8,
                                                      rhos = (0.7, 0.8, 0.9))
        # test conservion between primitive and conservative variables
        prim_vars = SVector(H..., v1..., b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
    end

    let equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 9.8,
                                                      rhos = (0.7, 0.8, 0.9))
        # test conservion between primitive and conservative variables
        prim_vars = SVector(H..., v1..., v2..., b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

        total_energy = energy_total(cons_vars, equations)
        @test total_energy ≈ entropy(cons_vars, equations)
        @test total_energy ≈
              energy_internal(cons_vars, equations) +
              energy_kinetic(cons_vars, equations)
    end
end

@timed_testset "SWE-Exner conversion between conservative/entropy variables" begin
    h, v, h_b = (1.0, 0.3, 0.1)

    let equations = ShallowWaterExnerEquations1D(gravity_constant = 9.81, rho_f = 0.9,
                                                 rho_s = 1.0, porosity = 0.4,
                                                 sediment_model = GrassModel(A_g = 0.01))
        # Test conversion between primitive and conservative variables
        prim_vars = SVector(h, v, h_b)
        cons_vars = prim2cons(prim_vars, equations)
        @test prim_vars ≈ cons2prim(cons_vars, equations)

        # Test conversion from conservative to entropy variables
        entropy_vars = cons2entropy(cons_vars, equations)
        @test entropy_vars ≈
              Trixi.ForwardDiff.gradient(u -> entropy(u, equations), cons_vars)
    end
end

@timed_testset "SWE-Exner eigenvalue / eigenvector computation" begin
    h, v, h_b = (1.0, 0.3, 0.1)

    let equations = ShallowWaterExnerEquations1D(gravity_constant = 9.81, rho_f = 0.9,
                                                 rho_s = 1.0, porosity = 0.4,
                                                 sediment_model = GrassModel(A_g = 0.01))
        u = SVector(h, h * v, h_b)
        r = equations.r
        g = equations.gravity

        # Compute effective sediment height
        h_s = TrixiShallowWater.q_s(SVector(h, h * v, 0.0), equations) / v
        dq_s_dh, dq_s_dhv, _ = Trixi.ForwardDiff.gradient(u -> TrixiShallowWater.q_s(u,
                                                                                     equations),
                                                          u)

        # flux Jacobian
        A = [0 1 0; (g * (h + h_s)-v^2) (2*v) (g*(h + 1 / equations.r * h_s));
             dq_s_dh dq_s_dhv 0]

        # Compute the eigenvalues using Cardano's formula
        λ1, λ2, λ3 = TrixiShallowWater.eigvals_cardano(SVector(h, h * v, h_b),
                                                       equations)

        # Precompute some common expressions
        c1 = g * (h + h_s)
        c2 = g * (h + h_s / r)

        # Eigenvector matrix
        r31 = ((v - λ1)^2 - c1) / c2
        r32 = ((v - λ2)^2 - c1) / c2
        r33 = ((v - λ3)^2 - c1) / c2
        R = [[1 1 1]; [λ1 λ2 λ3]; [r31 r32 r33]]

        # Inverse eigenvector matrix
        d1 = (λ1 - λ2) * (λ1 - λ3)
        d2 = (λ2 - λ1) * (λ2 - λ3)
        d3 = (λ3 - λ2) * (λ3 - λ1)
        R_inv = [(c1 - v^2 + λ2 * λ3)/d1 (2 * v - λ2 - λ3)/d1 c2/d1;
                 (c1 - v^2 + λ1 * λ3)/d2 (2 * v - λ1 - λ3)/d2 c2/d2;
                 (c1 - v^2 + λ1 * λ2)/d3 (2 * v - λ2 - λ1)/d3 c2/d3]

        # Eigenvalue vale matrix
        Λ = [λ1 0 0; 0 λ2 0; 0 0 λ3]

        @test R * R_inv ≈ [1 0 0; 0 1 0; 0 0 1]
        @test A ≈ R * Λ * R_inv
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
                                                                     0.2
                                                                 ])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81,
                                                                 rhos = [
                                                                     -1.0,
                                                                     0.1,
                                                                     0.2
                                                                 ])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                                                 rhos = [0.1, 0.3, 0.2])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81,
                                                                 rhos = [0.1, 0.3, 0.2])
    # Ensure that both tuple and array input are equivalent    
    @test ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                            rhos = [0.1, 0.2, 0.3]) ==
          ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                            rhos = (0.1, 0.2, 0.3))
    @test ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81,
                                            rhos = [0.1, 0.2, 0.3]) ==
          ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81,
                                            rhos = (0.1, 0.2, 0.3))

    # Check for argument error if initial condition is called with wrong number of layers
    equations = ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                                  rhos = (0.1, 0.2, 0.3, 0.4))
    @test_throws ArgumentError initial_condition_convergence_test(0.0, 0.0, equations)
    equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81,
                                                  rhos = (0.1, 0.2, 0.3, 0.4))
    @test_throws ArgumentError initial_condition_convergence_test(0.0, 0.0, equations)
end

@timed_testset "Exception check for default_threshold functions" begin
    @test_throws ArgumentError TrixiShallowWater.default_threshold_partially_wet(Int64)
    @test_throws ArgumentError TrixiShallowWater.default_threshold_desingularization(Int64)
end
end # Unit tests

end # module
