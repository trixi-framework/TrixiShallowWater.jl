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
    using Test: @test_nowarn
    # OBS! Constructing indicators/controllers using the parameters below doesn't make sense. It's
    # just useful to run basic tests of `show` methods.

    indicator_hg_swe = IndicatorHennemannGassnerShallowWater(1.0, 0.0, true, "variable",
                                                             "cache")
    @test_nowarn show(stdout, indicator_hg_swe)
end

@testset "Conversion between conservative / entropy variables" begin
    @timed_testset "ShallowWaterEquations" begin
        H, v1, v2, b = 3.5, 0.25, 0.1, 0.4

        let equations = ShallowWaterEquations1D(gravity = 9.8)
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

        let equations = ShallowWaterEquations2D(gravity = 9.8)
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

    @timed_testset "ShallowWaterEquationsQuasi1D" begin
        H, v1, v2, b, a = 3.5, 0.25, 0.1, 0.4, 0.3
        let equations = ShallowWaterEquationsQuasi1D(gravity = 9.8)
            cons_vars = prim2cons(SVector(H, v1, b, a), equations)
            entropy_vars = cons2entropy(cons_vars, equations)

            total_energy = energy_total(cons_vars, equations)
            @test entropy(cons_vars, equations) ≈ a * total_energy
        end
    end

    @timed_testset "ShallowWaterTwoLayerEquations" begin
        H_upper, v1_upper, v2_upper, H_lower, v1_lower, v2_lower, b = 3.5, 0.25, 0.13,
                                                                      2.5,
                                                                      0.1, 0.37, 0.4

        let equations = ShallowWaterTwoLayerEquations1D(gravity = 9.8,
                                                        rho_upper = 0.9,
                                                        rho_lower = 1.0)
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

        let equations = ShallowWaterTwoLayerEquations2D(gravity = 9.8,
                                                        rho_upper = 0.9,
                                                        rho_lower = 1.0)
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

    @timed_testset "ShallowWaterMultiLayerEquations" begin
        H = (3.5, 2.5, 1.5)
        v1 = (0.25, 0.1, 0.37)
        v2 = (0.13, 0.2, 0.3)
        b = 0.4

        let equations = ShallowWaterMultiLayerEquations1D(gravity = 9.8,
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

        let equations = ShallowWaterMultiLayerEquations2D(gravity = 9.8,
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

    @timed_testset "ShallowWaterExnerEquations" begin
        h, v, h_b = (1.0, 0.3, 0.1)

        let equations = ShallowWaterExnerEquations1D(gravity = 9.81, rho_f = 0.9,
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
end

@timed_testset "Consistency check for HLL flux (naive): SWE" begin
    flux_hll = FluxHLL(min_max_speed_naive)

    equations = ShallowWaterEquations1D(gravity = 9.81)
    u = SVector(1, 0.5, 0.0)
    @test flux_hll(u, u, 1, equations) ≈ flux(u, 1, equations)

    u_ll = SVector(0.1, 1.0, 0.0)
    u_rr = SVector(0.1, 1.0, 0.0)
    @test flux_hll(u_ll, u_rr, 1, equations) ≈ flux(u_ll, 1, equations)

    u_ll = SVector(0.1, -1.0, 0.0)
    u_rr = SVector(0.1, -1.0, 0.0)
    @test flux_hll(u_ll, u_rr, 1, equations) ≈ flux(u_rr, 1, equations)

    equations = ShallowWaterEquations2D(gravity = 9.81)
    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]
    u = SVector(1, 0.5, 0.5, 0.0)
    for normal_direction in normal_directions
        @test flux_hll(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    normal_direction = SVector(1.0, 0.0, 0.0)
    u_ll = SVector(0.1, 1.0, 1.0, 0.0)
    u_rr = SVector(0.1, 1.0, 1.0, 0.0)
    @test flux_hll(u_ll, u_rr, normal_direction, equations) ≈
          flux(u_ll, normal_direction, equations)

    u_ll = SVector(0.1, -1.0, -1.0, 0.0)
    u_rr = SVector(0.1, -1.0, -1.0, 0.0)
    @test flux_hll(u_ll, u_rr, normal_direction, equations) ≈
          flux(u_rr, normal_direction, equations)
end

@timed_testset "Consistency check for HLL flux with Davis wave speed estimates: SWE" begin
    flux_hll = FluxHLL(min_max_speed_davis)

    equations = ShallowWaterEquations1D(gravity = 9.81)
    u = SVector(1, 0.5, 0.0)
    @test flux_hll(u, u, 1, equations) ≈ flux(u, 1, equations)

    equations = ShallowWaterEquations2D(gravity = 9.81)
    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]
    u = SVector(1, 0.5, 0.5, 0.0)
    for normal_direction in normal_directions
        @test flux_hll(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    orientations = [1, 2]
    for orientation in orientations
        @test flux_hll(u, u, orientation, equations) ≈ flux(u, orientation, equations)
    end
end

@timed_testset "Consistency check for HLLE flux: SWE" begin
    equations = ShallowWaterEquations1D(gravity = 9.81)
    u = SVector(1, 0.5, 0.0)
    @test flux_hlle(u, u, 1, equations) ≈ flux(u, 1, equations)

    equations = ShallowWaterEquations2D(gravity = 9.81)
    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]
    orientations = [1, 2]

    u = SVector(1, 0.5, 0.5, 0.0)

    for orientation in orientations
        @test flux_hlle(u, u, orientation, equations) ≈ flux(u, orientation, equations)
    end

    for normal_direction in normal_directions
        @test flux_hlle(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end
end

@testset "Velocity functions for different equations" begin
    v1, v2 = pi, exp(1.0) # use pi, exp to test with non-trivial numbers
    v_vector = SVector(v1, v2)
    normal_direction_2d = SVector(pi^2, pi^3)
    v_normal_1d = v1 * normal_direction_2d[1]
    v_normal_2d = v1 * normal_direction_2d[1] + v2 * normal_direction_2d[2]
    H, b = exp(pi), exp(pi^2)
    gravity, H0 = 9.91, 0.1 # Standard numbers + 0.1

    shallow_water_1d = ShallowWaterEquations1D(; gravity, H0)
    u = prim2cons(SVector(H, v1, b), shallow_water_1d)
    @test isapprox(velocity(u, shallow_water_1d), v1)

    shallow_water_2d = ShallowWaterEquations2D(; gravity, H0)
    u = prim2cons(SVector(H, v1, v2, b), shallow_water_2d)
    @test isapprox(velocity(u, shallow_water_2d), SVector(v1, v2))
    @test isapprox(velocity(u, normal_direction_2d, shallow_water_2d), v_normal_2d)
    for orientation in 1:2
        @test isapprox(velocity(u, orientation, shallow_water_2d),
                       v_vector[orientation])
    end
end

@timed_testset "SWE-Exner eigenvalue / eigenvector computation" begin
    h, v, h_b = (1.0, 0.3, 0.1)

    let equations = ShallowWaterExnerEquations1D(gravity = 9.81, rho_f = 0.9,
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

    let equations = ShallowWaterEquations1D(gravity = 9.8)
        cons_vars = prim2cons(SVector(H, v1, b), equations)
        @test waterheight_pressure(cons_vars, equations) ≈
              waterheight(cons_vars, equations) * pressure(cons_vars, equations)
    end

    let equations = ShallowWaterEquations2D(gravity = 9.8)
        cons_vars = prim2cons(SVector(H, v1, v2, b), equations)
        @test waterheight_pressure(cons_vars, equations) ≈
              waterheight(cons_vars, equations) * pressure(cons_vars, equations)
    end
end

@testset "Equivalent Wave Speed Estimates: max_abs_speed(naive)" begin
    @timed_testset "ShallowWaterEquations1D" begin
        equations = ShallowWaterEquations1D(gravity = 9.81)

        h_ll_rr = SVector(12.0, 12.0)
        hv_ll_rr = SVector(42.0, 24.0)
        b_ll_rr = SVector(pi, pi)

        u_ll = SVector(h_ll_rr[1], hv_ll_rr[1], b_ll_rr[1])
        u_rr = SVector(h_ll_rr[2], hv_ll_rr[2], b_ll_rr[2])

        @test max_abs_speed_naive(u_ll, u_rr, 1, equations) ≈
              max_abs_speed(u_ll, u_rr, 1, equations)
    end

    @timed_testset "ShallowWaterEquations2D" begin
        equations = ShallowWaterEquations2D(gravity = 9.81)

        h_ll_rr = SVector(12.0, 12.0)
        hv1_ll_rr = SVector(42.0, 24.0)
        hv2_ll_rr = SVector(24.0, 42.0)
        b_ll_rr = SVector(pi, pi)

        u_ll = SVector(h_ll_rr[1], hv1_ll_rr[1], hv2_ll_rr[1], b_ll_rr[1])
        u_rr = SVector(h_ll_rr[2], hv1_ll_rr[2], hv2_ll_rr[2], b_ll_rr[2])

        for orientation in [1, 2]
            @test max_abs_speed_naive(u_ll, u_rr, orientation, equations) ≈
                  max_abs_speed(u_ll, u_rr, orientation, equations)
        end

        normal_directions = [SVector(1.0, 0.0),
            SVector(0.0, 1.0),
            SVector(0.5, -0.5),
            SVector(-1.2, 0.3)]

        for normal_direction in normal_directions
            @test max_abs_speed_naive(u_ll, u_rr, normal_direction, equations) ≈
                  max_abs_speed(u_ll, u_rr, normal_direction, equations)
        end
    end

    @timed_testset "ShallowWaterEquationsQuasi1D" begin
        equations = ShallowWaterEquationsQuasi1D(gravity = 9.81)

        ah_ll_rr = SVector(12.0, 12.0)
        ahv_ll_rr = SVector(42.0, 24.0)
        b_ll_rr = SVector(pi, pi)
        a_ll_rr = SVector(0.1, 0.1)

        u_ll = SVector(ah_ll_rr[1], ahv_ll_rr[1], b_ll_rr[1], a_ll_rr[1])
        u_rr = SVector(ah_ll_rr[2], ahv_ll_rr[2], b_ll_rr[2], a_ll_rr[2])

        @test max_abs_speed_naive(u_ll, u_rr, 1, equations) ≈
              max_abs_speed(u_ll, u_rr, 1, equations)
    end

    @timed_testset "ShallowWaterMultiLayerEquations1D" begin
        equations = ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                                      rhos = (1.0, 4.2))

        u_ll = SVector(1.0, 0.5, 0.6, 0.8)
        u_rr = SVector(1.0, 0.5, 0.3, 0.4)

        @test max_abs_speed_naive(u_ll, u_rr, 1, equations) ≈
              max_abs_speed(u_ll, u_rr, 1, equations)
    end

    @timed_testset "ShallowWaterMultiLayerEquations2D" begin
        equations = ShallowWaterMultiLayerEquations2D(gravity = 9.81,
                                                      rhos = (0.8, 1.2))

        u_ll = SVector(1.0, 1.0, 0.6, 0.6, 0.8, 0.8)
        u_rr = SVector(1.0, 1.0, 0.3, 0.3, 0.4, 0.4)

        for orientation in [1, 2]
            @test max_abs_speed_naive(u_ll, u_rr, orientation, equations) ≈
                  max_abs_speed(u_ll, u_rr, orientation, equations)
        end

        normal_directions = [SVector(1.0, 0.0),
            SVector(0.0, 1.0),
            SVector(0.5, -0.5),
            SVector(-1.2, 0.3)]

        for normal_direction in normal_directions
            @test max_abs_speed_naive(u_ll, u_rr, normal_direction, equations) ≈
                  max_abs_speed(u_ll, u_rr, normal_direction, equations)
        end
    end

    @timed_testset "ShallowWaterTwoLayerEquations1D" begin
        equations = ShallowWaterTwoLayerEquations1D(gravity = 9.8,
                                                    rho_upper = 0.9, rho_lower = 1.0)

        u_ll = SVector(1.0, 0.5, 0.6, 0.8, 42)
        u_rr = SVector(1.0, 0.3, 0.6, 0.4, 42)

        @test max_abs_speed_naive(u_ll, u_rr, 1, equations) ≈
              max_abs_speed(u_ll, u_rr, 1, equations)
    end

    @timed_testset "ShallowWaterTwoLayerEquations2D" begin
        equations = ShallowWaterTwoLayerEquations2D(gravity = 9.81,
                                                    rho_upper = 0.9, rho_lower = 1.0)

        u_ll = SVector(1.0, 1.0, 0.6, 0.6, 0.8, 0.8, 42)
        u_rr = SVector(1.0, 0.5, 0.3, 0.6, 0.4, 0.4, 42)

        for orientation in [1, 2]
            @test max_abs_speed_naive(u_ll, u_rr, orientation, equations) ≈
                  max_abs_speed(u_ll, u_rr, orientation, equations)
        end

        normal_directions = [SVector(1.0, 0.0),
            SVector(0.0, 1.0),
            SVector(0.5, -0.5),
            SVector(-1.2, 0.3)]

        for normal_direction in normal_directions
            @test max_abs_speed_naive(u_ll, u_rr, normal_direction, equations) ≈
                  max_abs_speed(u_ll, u_rr, normal_direction, equations)
        end
    end
end # Equivalent Wave Speed Estimates: max_abs_speed(naive)

@timed_testset "Exception check for the 2LSWE" begin
    error_message = "Invalid input: Densities must be chosen such that rho_upper < rho_lower"
    @test_throws error_message ShallowWaterTwoLayerEquations1D(gravity = 9.81,
                                                               rho_upper = 1.0,
                                                               rho_lower = 0.9)
    @test_throws error_message ShallowWaterTwoLayerEquations2D(gravity = 9.81,
                                                               rho_upper = 1.0,
                                                               rho_lower = 0.9)
end

@timed_testset "Input argument check for the ML-SWE" begin
    @test_throws ArgumentError ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                                                 rhos = [
                                                                     -1.0,
                                                                     0.1,
                                                                     0.2
                                                                 ])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations2D(gravity = 9.81,
                                                                 rhos = [
                                                                     -1.0,
                                                                     0.1,
                                                                     0.2
                                                                 ])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                                                 rhos = [0.1, 0.3, 0.2])
    @test_throws ArgumentError ShallowWaterMultiLayerEquations2D(gravity = 9.81,
                                                                 rhos = [0.1, 0.3, 0.2])
    # Ensure that both tuple and array input are equivalent
    @test ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                            rhos = [0.1, 0.2, 0.3]) ==
          ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                            rhos = (0.1, 0.2, 0.3))
    @test ShallowWaterMultiLayerEquations2D(gravity = 9.81,
                                            rhos = [0.1, 0.2, 0.3]) ==
          ShallowWaterMultiLayerEquations2D(gravity = 9.81,
                                            rhos = (0.1, 0.2, 0.3))

    # Check for argument error if initial condition is called with wrong number of layers
    equations = ShallowWaterMultiLayerEquations1D(gravity = 9.81,
                                                  rhos = (0.1, 0.2, 0.3, 0.4))
    @test_throws ArgumentError initial_condition_convergence_test(0.0, 0.0, equations)
    equations = ShallowWaterMultiLayerEquations2D(gravity = 9.81,
                                                  rhos = (0.1, 0.2, 0.3, 0.4))
    @test_throws ArgumentError initial_condition_convergence_test(0.0, 0.0, equations)
end

@timed_testset "Input argument check for SWE-Exner" begin
    # Type tests for GrassModel
    @test typeof(GrassModel(A_g = 0.1)) === typeof(GrassModel(A_g = 0.1, m_g = 3)) ===
          GrassModel{Float64}
    @test typeof(GrassModel(A_g = 0.1f0)) ===
          typeof(GrassModel(A_g = 0.1f0, m_g = 3)) === GrassModel{Float32}

    # Type tests for MeyerPeterMueller
    @test typeof(MeyerPeterMueller(theta_c = 0, d_s = 1e-3)) ===
          ShieldsStressModel{Float64}
    @test typeof(MeyerPeterMueller(theta_c = 0, d_s = 1.0f-3)) ===
          ShieldsStressModel{Float32}

    # Type tests for general ShieldsStressModel
    @test typeof(ShieldsStressModel(0.0, 1.5, 0.0, 8.0, 1.0, 0.0, 0.0, 1e-3)) ===
          ShieldsStressModel{Float64}
    @test typeof(ShieldsStressModel(0.0f0, 1.5f0, 0.0f0, 8.0f0, 1.0f0, 0.0f0, 0.0f0,
                                    1.0f-3)) === ShieldsStressModel{Float32}
    @test_throws MethodError ShieldsStressModel(0, 1.5, 0, 8, 1, 0, 0, 1e-3)
end

@timed_testset "Exception check for default_threshold functions" begin
    @test_throws ArgumentError TrixiShallowWater.default_threshold_partially_wet(Int64)
    @test_throws ArgumentError TrixiShallowWater.default_threshold_desingularization(Int64)
end

@timed_testset "Consistency check for boundary condition arguments" begin
    let
        equations = ShallowWaterEquations2D(gravity = 9.81)
        u_inner = SVector(1.0, 0.3, 0.3, 0.1)
        t = 1.0
        RealT = typeof(equations.gravity)

        # Check that the boundary condition functions have the correct output type
        boundary_condition = BoundaryConditionWaterHeight(1.0, equations)
        @test typeof(boundary_condition.h_boundary(t)) == RealT
        @test BoundaryConditionWaterHeight(t -> 1.0, equations).h_boundary(t) ==
              boundary_condition.h_boundary(t)
        @test BoundaryConditionWaterHeight(1, equations).h_boundary(t) ==
              boundary_condition.h_boundary(t)
        @test BoundaryConditionWaterHeight(1.0f0, equations).h_boundary(t) ==
              boundary_condition.h_boundary(t)

        @test_throws ArgumentError BoundaryConditionWaterHeight(t -> 1, equations)

        boundary_condition = BoundaryConditionMomentum(0.3, 0.1, equations)
        @test typeof(boundary_condition.hv_boundary(t)) == Tuple{RealT, RealT}
        @test BoundaryConditionMomentum(t -> 0.3, t -> 0.1, equations).hv_boundary(t) ==
              boundary_condition.hv_boundary(t)
        # Here we only check the type since 0.1f0 != 0.1
        @test typeof(BoundaryConditionMomentum(0.3, 0.1f0, equations).hv_boundary(t)) ==
              typeof(boundary_condition.hv_boundary(t))

        @test_throws ArgumentError BoundaryConditionMomentum(t -> 0.3 * t, t -> 1,
                                                             equations)
        @test_throws ArgumentError BoundaryConditionMomentum(t -> 0.3 * t, t -> 1.0f0,
                                                             equations)
    end

    let
        equations = ShallowWaterEquations1D(gravity = 9.81)
        u_inner = SVector(1.0, 0.3, 0.1)
        t = 1.0
        RealT = typeof(equations.gravity)

        # Check that the boundary condition functions have the correct output type
        boundary_condition = BoundaryConditionWaterHeight(1.0, equations)
        @test typeof(boundary_condition.h_boundary(t)) == RealT
        @test BoundaryConditionWaterHeight(t -> 1.0, equations).h_boundary(t) ==
              boundary_condition.h_boundary(t)
        @test BoundaryConditionWaterHeight(1, equations).h_boundary(t) ==
              boundary_condition.h_boundary(t)
        @test BoundaryConditionWaterHeight(1.0f0, equations).h_boundary(t) ==
              boundary_condition.h_boundary(t)

        @test_throws ArgumentError BoundaryConditionWaterHeight(t -> 1, equations)

        boundary_condition = BoundaryConditionMomentum(0.3, equations)
        @test typeof(boundary_condition.hv_boundary(t)) == RealT
        @test BoundaryConditionMomentum(t -> 0.3, equations).hv_boundary(t) ==
              boundary_condition.hv_boundary(t)
        # Here we only check the type since 0.1f0 != 0.1
        @test typeof(BoundaryConditionMomentum(0.3f0, equations).hv_boundary(t)) ==
              typeof(boundary_condition.hv_boundary(t))

        @test_throws ArgumentError BoundaryConditionMomentum(t -> 0.3f0,
                                                             equations)
    end
end
end # Unit tests

end # module
