module TestType

using Test
using Trixi
using TrixiShallowWater

include("test_trixi.jl")

# Start with a clean environment: remove Trixi.jl output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive = true)

@testset "Test Type Stability" begin
    @timed_testset "Shallow Water 1D" begin
        for RealT in (Float32, Float64)
            equations = @inferred ShallowWaterEquations1D(gravity = RealT(9.81),
                                                          threshold_partially_wet = RealT(1e-4),
                                                          threshold_desingularization = RealT(1e-6))

            x = SVector(zero(RealT))
            t = zero(RealT)
            u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT))
            orientation = 1
            directions = [1, 2]
            normal_direction = SVector(one(RealT))

            surface_flux_function = (flux_lax_friedrichs,
                                     flux_nonconservative_wintermeyer_etal)
            dissipation = DissipationLocalLaxFriedrichs()
            numflux = FluxHLL()

            @test eltype(@inferred initial_condition_convergence_test(x, t, equations)) ==
                  RealT
            @test eltype(@inferred initial_condition_weak_blast_wave(x, t, equations)) ==
                  RealT
            @test eltype(@inferred source_terms_convergence_test(u, x, t, equations)) ==
                  RealT

            for direction in directions
                @test eltype(@inferred boundary_condition_slip_wall(u_inner,
                                                                    orientation,
                                                                    direction,
                                                                    x, t,
                                                                    surface_flux_function,
                                                                    equations)) ==
                      SVector{3, RealT}
                @test eltype(@inferred TrixiShallowWater.calc_wavespeed_roe(u_ll, u_rr,
                                                                            direction,
                                                                            equations)) ==
                      RealT
            end

            @test eltype(@inferred flux(u, orientation, equations)) == RealT
            @test eltype(@inferred flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                                         orientation,
                                                                         equations)) ==
                  RealT
            @test eltype(@inferred flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                                       orientation,
                                                                       equations)) == RealT
            @test eltype(@inferred flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                                     orientation,
                                                                     equations)) == RealT
            @test eltype(@inferred flux_fjordholm_etal(u_ll, u_rr, orientation,
                                                       equations)) == RealT
            @test eltype(@inferred flux_wintermeyer_etal(u_ll, u_rr, orientation,
                                                         equations)) == RealT

            @test eltype(eltype(@inferred hydrostatic_reconstruction_audusse_etal(u_ll,
                                                                                  u_rr,
                                                                                  equations))) ==
                  RealT

            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation, equations)) ==
                  RealT
            @test eltype(@inferred dissipation(u_ll, u_rr, orientation, equations)) == RealT
            @test eltype(@inferred dissipation(u_ll, u_rr, normal_direction, equations)) ==
                  RealT
            @test eltype(@inferred numflux(u_ll, u_rr, orientation, equations)) == RealT
            # no matching method
            # @test eltype(@inferred numflux(u_ll, u_rr, normal_direction, equations)) == RealT
            @test eltype(@inferred min_max_speed_naive(u_ll, u_rr, orientation, equations)) ==
                  RealT
            @test eltype(@inferred min_max_speed_davis(u_ll, u_rr, orientation, equations)) ==
                  RealT
            @test eltype(@inferred min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                                          equations)) == RealT
            @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT

            @test typeof(@inferred velocity(u, equations)) == RealT
            @test eltype(@inferred cons2prim(u, equations)) == RealT
            @test eltype(@inferred prim2cons(u, equations)) == RealT
            @test eltype(@inferred cons2entropy(u, equations)) == RealT
            @test eltype(@inferred entropy2cons(u, equations)) == RealT
            @test typeof(@inferred waterheight(u, equations)) == RealT
            @test typeof(@inferred pressure(u, equations)) == RealT
            @test typeof(@inferred waterheight_pressure(u, equations)) == RealT

            @test typeof(@inferred entropy(cons, equations)) == RealT
            @test typeof(@inferred energy_total(cons, equations)) == RealT
            @test typeof(@inferred energy_kinetic(u, equations)) == RealT
            @test typeof(@inferred energy_internal(cons, equations)) == RealT
            # TODO: remove this?
            #@test typeof(@inferred lake_at_rest_error(u, equations)) == RealT
        end
    end

    @timed_testset "Shallow Water 2D" begin
        for RealT in (Float32, Float64)
            equations = @inferred ShallowWaterEquations2D(gravity = RealT(9.81),
                                                          threshold_partially_wet = RealT(1e-4),
                                                          threshold_desingularization = RealT(1e-6))

            x = SVector(zero(RealT), zero(RealT))
            t = zero(RealT)
            u = u_ll = u_rr = u_inner = cons = SVector(one(RealT), one(RealT), one(RealT),
                                                       one(RealT))
            orientations = [1, 2]
            directions = [1, 2, 3, 4]
            normal_direction = SVector(one(RealT), zero(RealT))

            surface_flux_function = (flux_lax_friedrichs,
                                     flux_nonconservative_wintermeyer_etal)
            dissipation = DissipationLocalLaxFriedrichs()
            numflux = FluxHLL()

            @test eltype(@inferred initial_condition_convergence_test(x, t, equations)) ==
                  RealT
            @test eltype(@inferred source_terms_convergence_test(u, x, t, equations)) ==
                  RealT
            @test eltype(@inferred initial_condition_weak_blast_wave(x, t, equations)) ==
                  RealT
            @test eltype(@inferred boundary_condition_slip_wall(u_inner,
                                                                normal_direction,
                                                                x, t,
                                                                surface_flux_function,
                                                                equations)) ==
                  SVector{4, RealT}

            @test eltype(@inferred velocity(u, normal_direction, equations)) == RealT
            @test eltype(@inferred flux(u, normal_direction, equations)) == RealT
            @test eltype(@inferred flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                                         normal_direction,
                                                                         equations)) ==
                  RealT
            @test eltype(@inferred flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                                       normal_direction,
                                                                       equations)) == RealT
            @test eltype(@inferred flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                                     normal_direction,
                                                                     equations)) == RealT
            @test eltype(@inferred flux_fjordholm_etal(u_ll, u_rr, normal_direction,
                                                       equations)) == RealT
            @test eltype(@inferred flux_wintermeyer_etal(u_ll, u_rr, normal_direction,
                                                         equations)) == RealT

            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, normal_direction,
                                                       equations)) ==
                  RealT
            @test eltype(@inferred dissipation(u_ll, u_rr, normal_direction, equations)) ==
                  RealT
            @test eltype(@inferred numflux(u_ll, u_rr, normal_direction, equations)) ==
                  RealT
            @test eltype(@inferred min_max_speed_naive(u_ll, u_rr, normal_direction,
                                                       equations)) == RealT
            @test eltype(@inferred min_max_speed_davis(u_ll, u_rr, normal_direction,
                                                       equations)) == RealT
            @test eltype(@inferred min_max_speed_einfeldt(u_ll, u_rr, normal_direction,
                                                          equations)) == RealT
            @test eltype(@inferred TrixiShallowWater.calc_wavespeed_roe(u_ll, u_rr,
                                                                        normal_direction,
                                                                        equations)) == RealT

            for orientation in orientations
                for direction in directions
                    @test eltype(@inferred boundary_condition_slip_wall(u_inner,
                                                                        orientation,
                                                                        direction, x, t,
                                                                        surface_flux_function,
                                                                        equations)) ==
                          SVector{4, RealT}
                end

                @test eltype(@inferred velocity(u, orientation, equations)) == RealT
                @test eltype(@inferred flux(u, orientation, equations)) == RealT
                @test eltype(@inferred flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                                             orientation,
                                                                             equations)) ==
                      RealT
                @test eltype(@inferred flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                                           orientation,
                                                                           equations)) ==
                      RealT
                @test eltype(@inferred flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                                         orientation,
                                                                         equations)) ==
                      RealT
                @test eltype(@inferred flux_fjordholm_etal(u_ll, u_rr, orientation,
                                                           equations)) == RealT
                @test eltype(@inferred flux_wintermeyer_etal(u_ll, u_rr, orientation,
                                                             equations)) == RealT

                @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation,
                                                           equations)) ==
                      RealT
                @test eltype(@inferred dissipation(u_ll, u_rr, orientation, equations)) ==
                      RealT
                @test eltype(@inferred numflux(u_ll, u_rr, orientation, equations)) == RealT
                @test eltype(@inferred min_max_speed_naive(u_ll, u_rr, orientation,
                                                           equations)) == RealT
                @test eltype(@inferred min_max_speed_davis(u_ll, u_rr, orientation,
                                                           equations)) == RealT
                @test eltype(@inferred min_max_speed_einfeldt(u_ll, u_rr, orientation,
                                                              equations)) == RealT
                @test eltype(@inferred TrixiShallowWater.calc_wavespeed_roe(u_ll, u_rr,
                                                                            orientation,
                                                                            equations)) ==
                      RealT
            end

            @test eltype(eltype(@inferred hydrostatic_reconstruction_audusse_etal(u_ll,
                                                                                  u_rr,
                                                                                  equations))) ==
                  RealT

            @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
            @test eltype(@inferred velocity(u, equations)) == RealT
            @test eltype(@inferred cons2prim(u, equations)) == RealT
            @test eltype(@inferred prim2cons(u, equations)) == RealT
            @test eltype(@inferred cons2entropy(u, equations)) == RealT
            @test eltype(@inferred entropy2cons(u, equations)) == RealT
            @test typeof(@inferred waterheight(u, equations)) == RealT
            @test typeof(@inferred pressure(u, equations)) == RealT
            @test typeof(@inferred waterheight_pressure(u, equations)) == RealT

            @test typeof(@inferred entropy(cons, equations)) == RealT
            @test typeof(@inferred energy_total(cons, equations)) == RealT
            @test typeof(@inferred energy_kinetic(u, equations)) == RealT
            @test typeof(@inferred energy_internal(cons, equations)) == RealT
            # TODO: remove this?
            #@test typeof(@inferred lake_at_rest_error(u, equations)) == RealT
        end
    end

    @timed_testset "Shallow Water Quasi 1D" begin
        for RealT in (Float32, Float64)
            equations = @inferred ShallowWaterEquationsQuasi1D(gravity = RealT(9.81))

            x = SVector(zero(RealT))
            t = zero(RealT)
            u = u_ll = u_rr = cons = SVector(one(RealT), one(RealT), one(RealT), one(RealT))
            orientation = 1
            normal_direction = normal_ll = normal_rr = SVector(one(RealT))

            dissipation = DissipationLocalLaxFriedrichs()

            @test eltype(@inferred initial_condition_convergence_test(x, t, equations)) ==
                  RealT
            @test eltype(@inferred source_terms_convergence_test(u, x, t, equations)) ==
                  RealT

            @test eltype(@inferred flux(u, orientation, equations)) == RealT
            @test eltype(@inferred flux_nonconservative_chan_etal(u_ll, u_rr,
                                                                  orientation,
                                                                  equations)) ==
                  RealT
            @test eltype(@inferred flux_nonconservative_chan_etal(u_ll, u_rr,
                                                                  normal_direction,
                                                                  equations)) ==
                  RealT
            @test eltype(@inferred flux_nonconservative_chan_etal(u_ll, u_rr, normal_ll,
                                                                  normal_rr,
                                                                  equations)) == RealT
            @test eltype(@inferred flux_chan_etal(u_ll, u_rr, orientation,
                                                  equations)) == RealT
            @test eltype(@inferred flux_chan_etal(u_ll, u_rr, normal_direction,
                                                  equations)) == RealT

            @test typeof(@inferred max_abs_speed_naive(u_ll, u_rr, orientation, equations)) ==
                  RealT
            @test eltype(@inferred dissipation(u_ll, u_rr, orientation, equations)) == RealT
            @test eltype(@inferred dissipation(u_ll, u_rr, normal_direction, equations)) ==
                  RealT
            @test eltype(@inferred Trixi.max_abs_speeds(u, equations)) == RealT
            @test typeof(@inferred velocity(u, equations)) == RealT
            @test eltype(@inferred cons2prim(u, equations)) == RealT
            @test eltype(@inferred prim2cons(u, equations)) == RealT
            @test eltype(@inferred cons2entropy(u, equations)) == RealT
            @test typeof(@inferred waterheight(u, equations)) == RealT

            @test typeof(@inferred entropy(cons, equations)) == RealT
            @test typeof(@inferred energy_total(cons, equations)) == RealT
            # TODO: remove this?
            #@test typeof(@inferred lake_at_rest_error(u, equations)) == RealT
        end
    end
end
end # module
