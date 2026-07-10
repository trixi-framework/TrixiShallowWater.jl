# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{2},
                                equations::ShallowWaterEquations2D, dg::DGSEM,
                                cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        # determine minimum value
        value_min = typemax(eltype(u))
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
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

        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)

            # Cut off velocity in case that the water height is smaller than the threshold.
            # Here the (possibly) cut off mean values are saved in a local variable
            # to ensure that it only influences the current node `i,j`.
            u_node, u_mean_local = zero_velocity_if_dry_node(u_node, u_mean, threshold,
                                                             equations)

            # When velocities are cut off, the only averaged value is the water height,
            # because the velocities are set to zero and this value is passed.
            # Otherwise, the velocities are averaged, as well.
            # Note that the auxiliary bottom topography variable `b` is never limited.
            set_node_vars!(u, theta * u_node + (1 - theta) * u_mean_local,
                           equations, dg, i, j, element)
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
                                           equations::ShallowWaterEquations2D)
    h_node, h_v1_node, h_v2_node, b_node = u_node
    # b_mean is not used as it must not be overwritten
    h_mean, h_v1_mean, h_v2_mean, _ = u_mean

    if h_node <= threshold
        h_v1_node = zero(eltype(u_node))
        h_v2_node = zero(eltype(u_node))
        h_v1_mean = zero(eltype(u_node))
        h_v2_mean = zero(eltype(u_node))
    end

    u_node = SVector(h_node, h_v1_node, h_v2_node, b_node)
    u_mean_local = SVector(h_mean, h_v1_mean, h_v2_mean, b_node)

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
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{2},
                                equations::ShallowWaterEquations2D, dg::DGSEM, cache,
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
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, new_element_id)
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
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, new_element_id)

                # Cut off velocity in case that the water height is smaller than the threshold.
                # Here the (possibly) cut off mean values are saved in a local variable
                # to ensure that it only influences the current node `i,j`.
                u_node, u_mean_local = zero_velocity_if_dry_node(u_node, u_mean,
                                                                 threshold, equations)

                # When velocities are cut off, the only averaged value is the water height,
                # because the velocities are set to zero and this value is passed.
                # Otherwise, the velocities are averaged, as well.
                # Note that the auxiliary bottom topography variable `b` is never limited.
                set_node_vars!(u, theta * u_node + (1 - theta) * u_mean_local,
                               equations, dg, i, j, new_element_id)
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
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{2},
                                equations::ShallowWaterEquations2D, dg::DGSEM, cache,
                                element_ids_new)
    # Apply limiter to coarsened elements
    Trixi.@threaded for element in element_ids_new
        # determine minimum value
        value_min = typemax(eltype(u))
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
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

        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)

            # Cut off velocity in case that the water height is smaller than the threshold.
            # Here the (possibly) cut off mean values are saved in a local variable
            # to ensure that it only influences the current node `i,j`.
            u_node, u_mean_local = zero_velocity_if_dry_node(u_node, u_mean, threshold,
                                                             equations)

            # When velocities are cut off, the only averaged value is the water height,
            # because the velocities are set to zero and this value is passed.
            # Otherwise, the velocities are averaged, as well.
            # Note that the auxiliary bottom topography variable `b` is never limited.
            set_node_vars!(u, theta * u_node + (1 - theta) * u_mean_local,
                           equations, dg, i, j, element)
        end
    end

    # "Safety" application of the wet/dry thresholds over all the DG nodes
    # on the current `element` after the limiting above in order to avoid dry nodes.
    # If the value_mean < threshold before applying limiter, there
    # could still be dry nodes afterwards due to logic of the limiting
    velocity_desingularization!(u, mesh, equations, dg, cache)

    return nothing
end

function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{2},
                                equations::ShallowWaterMultiLayerEquations2D, dg::DGSEM,
                                cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        # Limit layerwise
        for m in eachlayer(equations)
            # determine minimum value
            value_min = typemax(eltype(u))
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                value_min = min(value_min, variable(u_node, equations)[m])
            end

            # detect if limiting is necessary
            value_min < threshold - eps() || continue

            u_mean = Trixi.compute_u_mean(u, element, mesh, equations, dg, cache)

            # We compute the value directly with the mean values.
            # The waterheight `h` is limited independently in each layer.
            value_mean = variable(u_mean, equations)[m]
            theta = (value_mean - threshold) / (value_mean - value_min)

            # This avoids the issue when `value_mean` is slightly smaller than `threshold`
            # (e.g., due to finite precision effects in PositivityPreservingLimiterShallowWater),
            # which results in invalid theta values smaller than 0. Note that min(1, theta)
            # is not necessary since we are only enforcing lower bounds.
            theta = max(0, theta)

            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                h_node = waterheight(u_node, equations)[m]
                h_mean = waterheight(u_mean, equations)[m]

                u[m, i, j, element] = theta * h_node + (1 - theta) * h_mean
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

# Modified version of the limiter used in the refinement step of the AMR callback.
# To ensure admissibility after the refinement step, we compute a joint
# limiting coefficient for all children elements and then limit against the
# admissible mean value of the parent element.
# This strategy is described in Remark 3 of the paper:
# - Arpit Babbar, Praveen Chandrashekar (2025)
#   Lax-Wendroff flux reconstruction on adaptive curvilinear meshes with
#   error based time stepping for hyperbolic conservation laws
#   [doi: 10.1016/j.jcp.2024.113622](https://doi.org/10.1016/j.jcp.2024.113622)
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{2},
                                equations::ShallowWaterMultiLayerEquations2D, dg::DGSEM,
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
                for j in eachnode(dg), i in eachnode(dg)
                    u_node = get_node_vars(u, equations, dg, i, j, new_element_id)
                    value_min = min(value_min, variable(u_node, equations)[m])
                end
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
                for j in eachnode(dg), i in eachnode(dg)
                    u_node = get_node_vars(u, equations, dg, i, j, new_element_id)
                    h_node = waterheight(u_node, equations)[m]
                    h_mean = waterheight(u_mean, equations)[m]

                    u[m, i, j, new_element_id] = theta * h_node + (1 - theta) * h_mean
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
function limiter_shallow_water!(u, threshold::Real, variable,
                                mesh::Trixi.AbstractMesh{2},
                                equations::ShallowWaterMultiLayerEquations2D, dg::DGSEM,
                                cache,
                                element_ids_new)
    # Apply limiter to coarsened elements
    Trixi.@threaded for element in element_ids_new
        # Limit layerwise
        for m in eachlayer(equations)
            # determine minimum value
            value_min = typemax(eltype(u))
            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                value_min = min(value_min, variable(u_node, equations)[m])
            end
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

            for j in eachnode(dg), i in eachnode(dg)
                u_node = get_node_vars(u, equations, dg, i, j, element)
                h_node = waterheight(u_node, equations)[m]
                h_mean = waterheight(u_mean, equations)[m]

                u[m, i, j, element] = theta * h_node + (1 - theta) * h_mean
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
end # @muladd
