# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterEquations1D,
                                dg::DGSEM, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        # determine minimum value
        value_min = typemax(eltype(u))
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)
            value_min = min(value_min, variable(u_node, equations))
        end

        # detect if limiting is necessary
        value_min < threshold || continue

        u_mean = Trixi.compute_u_mean(u, element, mesh, equations, dg, cache)

        # We compute the value directly with the mean values, as we assume that
        # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
        # However, we only need linear variables such as the `waterheight` here,
        # where this is satisfied.
        value_mean = variable(u_mean, equations)
        theta = (value_mean - threshold) / (value_mean - value_min)

        # This avoids the issue when `value_mean` is slightly smaller than `threshold`
        # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
        # which results in invalid theta values smaller than 0. Note that min(1, theta)
        # is not necessary since we are only enforcing lower bounds.
        theta = max(0, theta)

        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)

            # Cut off velocity in case that the water height is smaller than the threshold.
            # Here the (possibly) cut off mean values are saved in a local variable
            # to ensure that it only influences the current node `i`.
            u_node, u_mean_local = zero_velocity_if_dry_node(u_node, u_mean, threshold,
                                                             equations)

            # When velocity is cut off, the only averaged value is the waterheight,
            # because the velocity is set to zero and this value is passed.
            # Otherwise, the velocity is averaged, as well.
            # Note that the auxiliary bottom topography variable `b` is never limited.
            set_node_vars!(u, theta * u_node + (1 - theta) * u_mean_local,
                           equations, dg, i, element)
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

# Cut off velocity (and thus the momenta) in case that the water height
# is smaller than the threshold.
@inline function zero_velocity_if_dry_node(u_node, u_mean, threshold,
                                           equations::ShallowWaterEquations1D)
    h_node, h_v_node, b_node = u_node
    # b_mean is not used as it must not be overwritten
    h_mean, h_v_mean, _ = u_mean

    if h_node <= threshold
        h_v_node = zero(eltype(u_node))
        h_v_mean = zero(eltype(u_node))
    end

    u_node = SVector(h_node, h_v_node, b_node)
    u_mean_local = SVector(h_mean, h_v_mean, b_node)

    return u_node, u_mean_local
end

# Modified version of the limiter used in the refinement step of the AMR callback.
# To ensure admissibility after the refinement step, we compute a joint
# limiting coefficient for all children elements and then limit against the
# admissible mean value of the parent element.
# This strategy is described in Remark 3 of the paper:
# - Arpit Babbar, Praveen Chandrashekar (2025)
#   Lax-Wendroff flux reconstruction on adaptive curvilinear meshes with
#   error based time stepping for hyperbolic conservation laws
#   [doi: 10.1016/j.jcp.2024.113622](https://doi.org/10.1016/j.jcp.2024.113622)
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterEquations1D, dg::DGSEM, cache,
                                element_ids_new, u_mean_refined_elements)
    @assert length(element_ids_new)==size(u_mean_refined_elements, 2) "The length of `element_ids_new` must match the second dimension of `u_mean_refined_elements`."

    Trixi.@threaded for idx in eachindex(element_ids_new)
        # Get the mean value from the parent element
        u_mean = get_node_vars(u_mean_refined_elements, equations, dg, idx)

        # We compute the value directly with the mean values, as we assume that
        # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
        # However, we only need linear variables such as the `waterheight` here,
        # where this is satisfied.
        value_mean = variable(u_mean, equations)
        theta = one(eltype(u)) # Limiting coefficient

        # Iterate over the children of the current element to determine a joint limiting coefficient `theta`
        for new_element_id in element_ids_new[idx]:(element_ids_new[idx] + 2^ndims(mesh) - 1)
            # determine minimum value
            value_min = typemax(eltype(u))
            for i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, new_element_id)
                value_min = min(value_min, variable(u_node, equations))
            end
            value_min < threshold || continue # Detect if limiting is necessary

            theta = min(theta, (value_mean - threshold) / (value_mean - value_min))
        end

        theta < 1 || continue # Check if limiting action is necessary

        # This avoids the issue when `value_mean` is slightly smaller than `threshold`
        # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
        # which results in invalid theta values smaller than 0. Note that min(1, theta)
        # is not necessary since we are only enforcing lower bounds.
        theta = max(0, theta)

        # Iterate again over the children to apply joint shifting
        for new_element_id in element_ids_new[idx]:(element_ids_new[idx] + 2^ndims(mesh) - 1)
            for i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, new_element_id)

                # Cut off velocity in case that the water height is smaller than the threshold.
                # Here the (possibly) cut off mean values are saved in a local variable
                # to ensure that it only influences the current node `i`.
                u_node, u_mean_local = zero_velocity_if_dry_node(u_node, u_mean,
                                                                 threshold, equations)

                # When velocities are cut off, the only averaged value is the water height,
                # because the velocities are set to zero and this value is passed.
                # Otherwise, the velocities are averaged, as well.
                # Note that the auxiliary bottom topography variable `b` is never limited.
                set_node_vars!(u, theta * u_node + (1 - theta) * u_mean_local,
                               equations, dg, i, new_element_id)
            end
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

# Modified version of the limiter used in the coarsening step of the AMR callback.
# To ensure admissibility after the coarsening step, we apply the limiter to
# the coarsened elements.
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterEquations1D, dg::DGSEM, cache,
                                element_ids_new)
    # Apply limiter to coarsened elements
    Trixi.@threaded for element in element_ids_new
        # determine minimum value
        value_min = typemax(eltype(u))
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)
            value_min = min(value_min, variable(u_node, equations))
        end
        value_min < threshold || continue # Detect if limiting is necessary

        u_mean = Trixi.compute_u_mean(u, element, mesh, equations, dg, cache)

        # We compute the value directly with the mean values, as we assume that
        # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
        # However, we only need linear variables such as the `waterheight` here,
        # where this is satisfied.
        value_mean = variable(u_mean, equations)
        theta = (value_mean - threshold) / (value_mean - value_min)

        # This avoids the issue when `value_mean` is slightly smaller than `threshold`
        # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
        # which results in invalid theta values smaller than 0. Note that min(1, theta)
        # is not necessary since we are only enforcing lower bounds.
        theta = max(0, theta)

        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)

            # Cut off velocity in case that the water height is smaller than the threshold.
            # Here the (possibly) cut off mean values are saved in a local variable
            # to ensure that it only influences the current node `i`.
            u_node, u_mean_local = zero_velocity_if_dry_node(u_node, u_mean, threshold,
                                                             equations)

            # When velocities are cut off, the only averaged value is the water height,
            # because the velocities are set to zero and this value is passed.
            # Otherwise, the velocities are averaged, as well.
            # Note that the auxiliary bottom topography variable `b` is never limited.
            set_node_vars!(u, theta * u_node + (1 - theta) * u_mean_local,
                           equations, dg, i, element)
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

# Note that for the `ShallowWaterMultiLayerEquations1D` only the waterheight `h` is limited in
# each layer. Furthermore, a velocity desingularization is applied after the limiting to avoid
# numerical problems near dry states.
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterMultiLayerEquations1D,
                                dg::DGSEM, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        # Limit layerwise
        for m in eachlayer(equations)
            # determine minimum value
            value_min = typemax(eltype(u))
            for i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, element)
                value_min = min(value_min, variable(u_node, equations)[m])
            end

            # detect if limiting is necessary
            # TODO: Determine default thresholds for Float32
            # (Currently default thresholds are only provided for Float64)
            value_min < threshold - eps() || continue

            u_mean = Trixi.compute_u_mean(u, element, mesh, equations, dg, cache)

            # We compute the value directly with the mean values.
            # Note: The waterheight `h` is limited independently in each layer.
            value_mean = variable(u_mean, equations)[m]
            theta = (value_mean - threshold) / (value_mean - value_min)

            # This avoids the issue when `value_mean` is slightly smaller than `threshold`
            # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
            # which results in invalid theta values smaller than 0. Note that min(1, theta)
            # is not necessary since we are only enforcing lower bounds.
            theta = max(0, theta)

            for i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, element)
                h_node = waterheight(u_node, equations)[m]
                h_mean = waterheight(u_mean, equations)[m]

                u[m, i, element] = theta * h_node + (1 - theta) * h_mean
            end
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting.
    # Additionally, a velocity desingularization is applied to avoid numerical
    # problems near dry states.
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

# Modified version of the limiter used in the refinement step of the AMR callback.
# To ensure admissibility after the refinement step, we compute a joint
# limiting coefficient for all children elements and then limit against the
# admissible mean value of the parent element.
# This strategy is described in Remark 3 of the paper:
# - Arpit Babbar, Praveen Chandrashekar (2025)
#   Lax-Wendroff flux reconstruction on adaptive curvilinear meshes with
#   error based time stepping for hyperbolic conservation laws
#   [doi: 10.1016/j.jcp.2024.113622](https://doi.org/10.1016/j.jcp.2024.113622)
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterMultiLayerEquations1D, dg::DGSEM,
                                cache,
                                element_ids_new, u_mean_refined_elements)
    @assert length(element_ids_new)==size(u_mean_refined_elements, 2) "The length of `element_ids_new` must match the second dimension of `u_mean_refined_elements`."

    Trixi.@threaded for idx in eachindex(element_ids_new)
        # Limit layerwise
        for m in eachlayer(equations)
            # Get the mean value from the parent element
            u_mean = get_node_vars(u_mean_refined_elements, equations, dg, idx)

            # We compute the value directly with the mean values, as we assume that
            # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
            # However, we only need linear variables such as the `waterheight` here,
            # where this is satisfied.
            # Note: The waterheight `h` is limited independently in each layer.
            value_mean = variable(u_mean, equations)[m]
            theta = one(eltype(u)) # Limiting coefficient

            # Iterate over the children of the current element to determine a joint limiting coefficient `theta`
            for new_element_id in element_ids_new[idx]:(element_ids_new[idx] + 2^ndims(mesh) - 1)
                # determine minimum value
                value_min = typemax(eltype(u))
                for i in eachnode(dg)
                    u_node = get_node_vars(u, equations, dg, i, new_element_id)
                    value_min = min(value_min, variable(u_node, equations)[m])
                end
                # TODO: Determine default thresholds for Float32
                # (Currently default thresholds are only provided for Float64)
                value_min < threshold - eps() || continue # Detect if limiting is necessary

                theta = min(theta, (value_mean - threshold) / (value_mean - value_min))
            end

            theta < 1 || continue # Check if limiting action is necessary

            # This avoids the issue when `value_mean` is slightly smaller than `threshold`
            # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
            # which results in invalid theta values smaller than 0. Note that min(1, theta)
            # is not necessary since we are only enforcing lower bounds.
            theta = max(0, theta)

            # Iterate again over the children to apply joint shifting
            for new_element_id in element_ids_new[idx]:(element_ids_new[idx] + 2^ndims(mesh) - 1)
                for i in eachnode(dg)
                    u_node = get_node_vars(u, equations, dg, i, new_element_id)
                    h_node = waterheight(u_node, equations)[m]
                    h_mean = waterheight(u_mean, equations)[m]

                    u[m, i, new_element_id] = theta * h_node + (1 - theta) * h_mean
                end
            end
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

# Modified version of the limiter used in the coarsening step of the AMR callback.
# To ensure admissibility after the coarsening step, we apply the limiter to
# the coarsened elements.
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterMultiLayerEquations1D, dg::DGSEM,
                                cache,
                                element_ids_new)
    # Apply limiter to coarsened elements
    Trixi.@threaded for element in element_ids_new
        # Limit layerwise
        for m in eachlayer(equations)
            # determine minimum value
            value_min = typemax(eltype(u))
            for i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, element)
                value_min = min(value_min, variable(u_node, equations)[m])
            end
            # TODO: Determine default thresholds for Float32
            # (Currently default thresholds are only provided for Float64)
            value_min < threshold - eps() || continue # Detect if limiting is necessary

            u_mean = Trixi.compute_u_mean(u, element, mesh, equations, dg, cache)

            # We compute the value directly with the mean values, as we assume that
            # Jensen's inequality holds (e.g. pressure for compressible Euler equations).
            # However, we only need linear variables such as the `waterheight` here,
            # where this is satisfied.
            # Note: The waterheight `h` is limited independently in each layer.
            value_mean = variable(u_mean, equations)[m]
            theta = (value_mean - threshold) / (value_mean - value_min)

            # This avoids the issue when `value_mean` is slightly smaller than `threshold`
            # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
            # which results in invalid theta values smaller than 0. Note that min(1, theta)
            # is not necessary since we are only enforcing lower bounds.
            theta = max(0, theta)

            for i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, element)
                h_node = waterheight(u_node, equations)[m]
                h_mean = waterheight(u_mean, equations)[m]

                u[m, i, element] = theta * h_node + (1 - theta) * h_mean
            end
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{1},
                                equations::ShallowWaterExnerEquations1D,
                                dg::DGSEM, cache)
    @unpack weights = dg.basis

    Trixi.@threaded for element in eachelement(dg, cache)
        # determine minimum value
        value_min = typemax(eltype(u))
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)
            value_min = min(value_min, variable(u_node, equations))
        end

        # detect if limiting is necessary
        value_min < threshold - eps() || continue

        # compute mean value
        u_mean = zero(get_node_vars(u, equations, dg, 1, element))
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)
            u_mean += u_node * weights[i]
        end

        # note that the reference element is [-1,1]^ndims(dg), thus the weights sum to 2
        u_mean = u_mean / 2^ndims(mesh)

        # We compute the value directly with the mean values.
        value_mean = variable(u_mean, equations)
        theta = (value_mean - threshold) / (value_mean - value_min)

        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)
            h_node = waterheight(u_node, equations)
            h_mean = waterheight(u_mean, equations)

            u[1, i, element] = theta * h_node + (1 - theta) * h_mean
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting.
    # Additionally, a velocity desingularization is applied to avoid numerical 
    # problems near dry states.
    Trixi.@threaded for element in eachelement(dg, cache)
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)
            h, hv, h_b = u_node

            # Velocity desingularization
            hv = h * (2.0 * h * hv) /
                 (h^2 + max(h^2, equations.threshold_desingularization))

            # Ensure positivity and zero velocity at dry states
            if h <= threshold
                h = threshold
                hv = zero(eltype(u))
            end

            u_node = SVector(h, hv, h_b)

            set_node_vars!(u, u_node, equations, dg, i, element)
        end
    end

    return nothing
end
end # @muladd
