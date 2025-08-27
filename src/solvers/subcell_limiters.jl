
# Specialized subcell limiter for shallow water equations.
# This version of the subcell limiter includes a modified alpha calculation to treat wet/dry elements.
# The idea is to set the limiter to pure FV in elements with a water height below a certain threshold.
# This is done by setting the limiter coefficient `alpha` to 1 in these elements.
# The threshold value can be set via `threshold_partially_wet` in the equations struct.
function (limiter::Trixi.SubcellLimiterIDP)(u::AbstractArray{<:Any, 4},
                                            semi,
                                            equations::Union{Trixi.AbstractShallowWaterEquations,
                                                             AbstractShallowWaterMultiLayerEquations},
                                            dg::DGSEM,
                                            t, dt;
                                            kwargs...)
    @unpack alpha = limiter.cache.subcell_limiter_coefficients
    # TODO: Do not abuse `reset_du!` but maybe implement a generic `set_zero!`
    Trixi.@trixi_timeit Trixi.timer() "reset alpha" Trixi.reset_du!(alpha, dg, semi.cache)

    if limiter.local_twosided
        Trixi.@trixi_timeit Trixi.timer() "local twosided" Trixi.idp_local_twosided!(alpha,
                                                                                     limiter,
                                                                                     u, t,
                                                                                     dt,
                                                                                     semi)
    end
    if limiter.positivity
        Trixi.@trixi_timeit Trixi.timer() "positivity" Trixi.idp_positivity!(alpha, limiter,
                                                                             u, dt, semi)
    end
    if limiter.local_onesided
        Trixi.@trixi_timeit Trixi.timer() "local onesided" Trixi.idp_local_onesided!(alpha,
                                                                                     limiter,
                                                                                     u, t,
                                                                                     dt,
                                                                                     semi)
    end

    # Calculate alpha1 and alpha2
    @unpack alpha1, alpha2 = limiter.cache.subcell_limiter_coefficients
    Trixi.@threaded for element in eachelement(dg, semi.cache)
        for j in eachnode(dg), i in 2:nnodes(dg)
            alpha1[i, j, element] = max(alpha[i - 1, j, element], alpha[i, j, element])
        end
        for j in 2:nnodes(dg), i in eachnode(dg)
            alpha2[i, j, element] = max(alpha[i, j - 1, element], alpha[i, j, element])
        end
        alpha1[1, :, element] .= zero(eltype(alpha1))
        alpha1[nnodes(dg) + 1, :, element] .= zero(eltype(alpha1))
        alpha2[:, 1, element] .= zero(eltype(alpha2))
        alpha2[:, nnodes(dg) + 1, element] .= zero(eltype(alpha2))
    end

    # Modification for wet/dry elements
    Trixi.@threaded for element in eachelement(dg, semi.cache)

        # (Re-)set dummy variable for alpha_dry
        indicator_wet = 1

        for j in eachnode(dg), i in 1:nnodes(dg)
            h = waterheight(u[:, i, j, element], equations)

            # Set indicator to FV if water height is below the threshold
            if minimum(h) <= equations.threshold_partially_wet
                indicator_wet = 0
            end
        end

        if indicator_wet == 0   # element is dry
            alpha[:, :, element] .= one(eltype(alpha))
            alpha1[:, :, element] .= one(eltype(alpha1))
            alpha2[:, :, element] .= one(eltype(alpha2))
        end
        # Reset the magic edges
        alpha1[1, :, element] .= zero(eltype(alpha1))
        alpha1[nnodes(dg) + 1, :, element] .= zero(eltype(alpha1))
        alpha2[:, 1, element] .= zero(eltype(alpha2))
        alpha2[:, nnodes(dg) + 1, element] .= zero(eltype(alpha2))
    end
    return nothing
end
