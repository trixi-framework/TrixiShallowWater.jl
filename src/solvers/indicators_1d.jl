# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Modified indicator for ShallowWaterEquationsWetDry1D and ShallowWaterMultiLayerEquations1D to 
# apply full FV method on elements containing some "dry" LGL nodes. That is, if an element is 
# partially "wet" then it becomes a full FV element.
function (indicator_hg::IndicatorHennemannGassnerShallowWater)(u::AbstractArray{<:Any,
                                                                                3},
                                                               mesh,
                                                               equations::Union{ShallowWaterEquationsWetDry1D,
                                                                                ShallowWaterMultiLayerEquations1D},
                                                               dg::DGSEM, cache;
                                                               kwargs...)
    @unpack alpha_max, alpha_min, alpha_smooth, variable = indicator_hg
    @unpack alpha, alpha_tmp, indicator_threaded, modal_threaded = indicator_hg.cache
    # TODO: Taal refactor, when to `resize!` stuff changed possibly by AMR?
    #       Shall we implement `resize!(semi::AbstractSemidiscretization, new_size)`
    #       or just `resize!` whenever we call the relevant methods as we do now?
    resize!(alpha, nelements(dg, cache))
    if alpha_smooth
        resize!(alpha_tmp, nelements(dg, cache))
    end

    # magic parameters
    threshold = 0.5 * 10^(-1.8 * (nnodes(dg))^0.25)
    parameter_s = log((1 - 0.0001) / 0.0001)

    #=
    If the water height `h` at one LGL node is lower than `threshold_partially_wet`
    the indicator sets the element-wise blending factor alpha[element] = 1
    via the local variable `indicator_wet`. In turn, this ensures that a pure
    FV method is used in partially wet elements and guarantees the well-balanced property.

    Hard-coded cut-off value of `threshold_partially_wet = 1e-4` was determined through many numerical experiments.
    Overall idea is to increase robustness when computing the velocity on (nearly) dry elements which
    could be "dangerous" due to division of conservative variables, e.g., v = hv / h.
    Here, the impact of the threshold on the number of elements being updated with FV is not that
    significant. However, its impact on the robustness is very significant.
    The value can be seen as a trade-off between accuracy and stability.
    Well-balancedness of the scheme on partially wet elements with hydrostatic reconstruction
    can only be proven for the FV method (see Chen and Noelle).
    Therefore we set alpha to one regardless of its given maximum value. 
    =#
    threshold_partially_wet = equations.threshold_partially_wet

    Trixi.@threaded for element in eachelement(dg, cache)
        indicator = indicator_threaded[Threads.threadid()]
        modal = modal_threaded[Threads.threadid()]

        # (Re-)set dummy variable for alpha_dry
        indicator_wet = 1

        # Calculate indicator variables at Gauss-Lobatto nodes
        for i in eachnode(dg)
            u_local = get_node_vars(u, equations, dg, i, element)
            h = waterheight(u_local, equations)

            # Set indicator to FV if water height is below the threshold
            if minimum(h) <= threshold_partially_wet
                indicator_wet = 0
            end

            indicator[i] = indicator_hg.variable(u_local, equations)
        end

        # Convert to modal representation
        Trixi.multiply_scalar_dimensionwise!(modal,
                                             dg.basis.inverse_vandermonde_legendre,
                                             indicator)

        # Calculate total energies for all modes, without highest, without two highest
        total_energy = zero(eltype(modal))
        for i in 1:nnodes(dg)
            total_energy += modal[i]^2
        end
        total_energy_clip1 = zero(eltype(modal))
        for i in 1:(nnodes(dg) - 1)
            total_energy_clip1 += modal[i]^2
        end
        total_energy_clip2 = zero(eltype(modal))
        for i in 1:(nnodes(dg) - 2)
            total_energy_clip2 += modal[i]^2
        end

        # Calculate energy in higher modes
        energy = max((total_energy - total_energy_clip1) / total_energy,
                     (total_energy_clip1 - total_energy_clip2) / total_energy_clip1)

        alpha_element = 1 / (1 + exp(-parameter_s / threshold * (energy - threshold)))

        # Take care of the case close to pure DG
        if alpha_element < alpha_min
            alpha_element = zero(alpha_element)
        end

        # Take care of the case close to pure FV
        if alpha_element > 1 - alpha_min
            alpha_element = one(alpha_element)
        end

        # Clip the maximum amount of FV allowed or set to one depending on indicator_wet
        if indicator_wet == 0
            alpha[element] = 1
        else # Element is not defined as dry but wet
            alpha[element] = min(alpha_max, alpha_element)
        end
    end

    if alpha_smooth
        Trixi.apply_smoothing!(mesh, alpha, alpha_tmp, dg, cache)
    end

    return alpha
end
end # @muladd
