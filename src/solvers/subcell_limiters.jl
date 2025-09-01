
# Specialized subcell limiter for shallow water equations.
# This version of the subcell limiter includes a modified alpha calculation to treat wet/dry elements.
# The idea is to set the limiter to pure FV in elements with a water height below a certain threshold.
# This is done by setting the limiter coefficient `alpha` to 1 in these elements.
# The threshold value can be set via `threshold_partially_wet` in the equations struct.
function (limiter::Trixi.SubcellLimiterIDP)(u::AbstractArray{<:Any, 4},
                                            semi,
                                            equations::AbstractShallowWaterMultiLayerEquations,
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

###############################################################################
# Calculation of local bounds using low-order FV solution.

# Specialized version of the two-sided limiter for the shallow water equations that uses the total 
# water height `H = h + b` for the limiting variable instead of the water height `h`.
# TODO: Add support for other mesh types. Right now only TreeMesh2D and P4estMesh2D are supported.
@inline function Trixi.calc_bounds_twosided!(var_min, var_max, variable,
                                             u, t, semi,
                                             equations::AbstractShallowWaterMultiLayerEquations{2,
                                                                                                4,
                                                                                                1})
    mesh, _, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    # Calc bounds inside elements
    Trixi.@threaded for element in eachelement(dg, cache)
        var_min[:, :, element] .= typemax(eltype(var_min))
        var_max[:, :, element] .= typemin(eltype(var_max))
        # Calculate bounds at Gauss-Lobatto nodes using u
        for j in eachnode(dg), i in eachnode(dg)
            var = u[variable, i, j, element]

            # If the water height is specified as the limiting variable, apply limiting with the 
            # total water height `H = h + b` instead.
            if variable == 1
                var += u[end, i, j, element]
            end

            var_min[i, j, element] = min(var_min[i, j, element], var)
            var_max[i, j, element] = max(var_max[i, j, element], var)

            if i > 1
                var_min[i - 1, j, element] = min(var_min[i - 1, j, element], var)
                var_max[i - 1, j, element] = max(var_max[i - 1, j, element], var)
            end
            if i < nnodes(dg)
                var_min[i + 1, j, element] = min(var_min[i + 1, j, element], var)
                var_max[i + 1, j, element] = max(var_max[i + 1, j, element], var)
            end
            if j > 1
                var_min[i, j - 1, element] = min(var_min[i, j - 1, element], var)
                var_max[i, j - 1, element] = max(var_max[i, j - 1, element], var)
            end
            if j < nnodes(dg)
                var_min[i, j + 1, element] = min(var_min[i, j + 1, element], var)
                var_max[i, j + 1, element] = max(var_max[i, j + 1, element], var)
            end
        end
    end

    # Values at element boundary
    Trixi.calc_bounds_twosided_interface!(var_min, var_max, variable,
                                          u, t, semi, mesh, equations)

    # If the water height is specified as the limiting variable, subtract the bottom topography
    # to get back the water height `h = H - b`
    if variable == 1
        Trixi.@threaded for element in eachelement(dg, cache)
            for j in eachnode(dg), i in eachnode(dg)
                var_min[i, j, element] -= u[end, i, j, element]
                var_max[i, j, element] -= u[end, i, j, element]
            end
        end
    end

    return nothing
end

@inline function Trixi.calc_bounds_twosided_interface!(var_min, var_max, variable,
                                                       u, t, semi, mesh::Trixi.TreeMesh2D,
                                                       equations::AbstractShallowWaterMultiLayerEquations{2,
                                                                                                          4,
                                                                                                          1})
    _, _, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    (; boundary_conditions) = semi
    # Calc bounds at interfaces and periodic boundaries
    for interface in Trixi.eachinterface(dg, cache)
        # Get neighboring element ids
        left = cache.interfaces.neighbor_ids[1, interface]
        right = cache.interfaces.neighbor_ids[2, interface]

        orientation = cache.interfaces.orientations[interface]

        for i in eachnode(dg)
            index_left = (nnodes(dg), i)
            index_right = (1, i)
            if orientation == 2
                index_left = reverse(index_left)
                index_right = reverse(index_right)
            end

            var_left = u[variable, index_left..., left]
            var_right = u[variable, index_right..., right]

            # If the water height is specified as the limiting variable, apply limiting with the 
            # total water height `H = h + b` instead.
            if variable == 1
                var_left += u[end, index_left..., left]
                var_right += u[end, index_right..., right]
            end

            var_min[index_right..., right] = min(var_min[index_right..., right],
                                                 var_left)
            var_max[index_right..., right] = max(var_max[index_right..., right],
                                                 var_left)

            var_min[index_left..., left] = min(var_min[index_left..., left], var_right)
            var_max[index_left..., left] = max(var_max[index_left..., left], var_right)
        end
    end

    # Calc bounds at physical boundaries
    for boundary in Trixi.eachboundary(dg, cache)
        element = cache.boundaries.neighbor_ids[boundary]
        orientation = cache.boundaries.orientations[boundary]
        neighbor_side = cache.boundaries.neighbor_sides[boundary]

        for i in eachnode(dg)
            if neighbor_side == 2 # Element is on the right, boundary on the left
                index = (1, i)
                boundary_index = 1
            else # Element is on the left, boundary on the right
                index = (nnodes(dg), i)
                boundary_index = 2
            end
            if orientation == 2
                index = reverse(index)
                boundary_index += 2
            end
            u_inner = get_node_vars(u, equations, dg, index..., element)
            u_outer = Trixi.get_boundary_outer_state(u_inner, t,
                                                     boundary_conditions[boundary_index],
                                                     orientation, boundary_index,
                                                     mesh, equations, dg, cache,
                                                     index..., element)
            var_outer = u_outer[variable]
            # If the water height is specified as the limiting variable, apply limiting with the
            # total water height `H = h + b` instead.
            if variable == 1
                var_outer += u_outer[end]
            end

            var_min[index..., element] = min(var_min[index..., element], var_outer)
            var_max[index..., element] = max(var_max[index..., element], var_outer)
        end
    end

    return nothing
end

function Trixi.calc_bounds_twosided_interface!(var_min, var_max, variable, u, t, semi,
                                               mesh::Trixi.P4estMesh{2},
                                               equations::AbstractShallowWaterMultiLayerEquations{2,
                                                                                                  4,
                                                                                                  1})
    _, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
    (; boundary_conditions) = semi

    (; neighbor_ids, node_indices) = cache.interfaces
    index_range = eachnode(dg)

    # Calc bounds at interfaces and periodic boundaries
    for interface in Trixi.eachinterface(dg, cache)
        # Get element and side index information on the primary element
        primary_element = neighbor_ids[1, interface]
        primary_indices = node_indices[1, interface]

        # Get element and side index information on the secondary element
        secondary_element = neighbor_ids[2, interface]
        secondary_indices = node_indices[2, interface]

        # Create the local i,j indexing
        i_primary_start, i_primary_step = Trixi.index_to_start_step_2d(primary_indices[1],
                                                                       index_range)
        j_primary_start, j_primary_step = Trixi.index_to_start_step_2d(primary_indices[2],
                                                                       index_range)
        i_secondary_start, i_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[1],
                                                                           index_range)
        j_secondary_start, j_secondary_step = Trixi.index_to_start_step_2d(secondary_indices[2],
                                                                           index_range)

        i_primary = i_primary_start
        j_primary = j_primary_start
        i_secondary = i_secondary_start
        j_secondary = j_secondary_start

        for node in eachnode(dg)
            var_primary = u[variable, i_primary, j_primary, primary_element]
            var_secondary = u[variable, i_secondary, j_secondary, secondary_element]

            # If the water height is specified as the limiting variable, apply limiting with the
            # total water height `H = h + b` instead.
            var_primary += u[end, i_primary, j_primary, primary_element]
            var_secondary += u[end, i_secondary, j_secondary, secondary_element]

            var_min[i_primary, j_primary, primary_element] = min(var_min[i_primary,
                                                                         j_primary,
                                                                         primary_element],
                                                                 var_secondary)
            var_max[i_primary, j_primary, primary_element] = max(var_max[i_primary,
                                                                         j_primary,
                                                                         primary_element],
                                                                 var_secondary)

            var_min[i_secondary, j_secondary, secondary_element] = min(var_min[i_secondary,
                                                                               j_secondary,
                                                                               secondary_element],
                                                                       var_primary)
            var_max[i_secondary, j_secondary, secondary_element] = max(var_max[i_secondary,
                                                                               j_secondary,
                                                                               secondary_element],
                                                                       var_primary)

            # Increment primary element indices
            i_primary += i_primary_step
            j_primary += j_primary_step
            i_secondary += i_secondary_step
            j_secondary += j_secondary_step
        end
    end

    # Calc bounds at physical boundaries
    Trixi.calc_bounds_twosided_boundary!(var_min, var_max, variable, u, t,
                                         boundary_conditions,
                                         mesh, equations, dg, cache)

    return nothing
end

@inline function Trixi.calc_bounds_twosided_boundary!(var_min, var_max, variable, u, t,
                                                      boundary_conditions::Trixi.BoundaryConditionPeriodic,
                                                      mesh::Trixi.P4estMesh{2},
                                                      equations::AbstractShallowWaterMultiLayerEquations{2,
                                                                                                         4,
                                                                                                         1},
                                                      dg, cache)
    return nothing
end

@inline function Trixi.calc_bounds_twosided_boundary!(var_min, var_max, variable, u, t,
                                                      boundary_conditions,
                                                      mesh::Trixi.P4estMesh{2},
                                                      equations::AbstractShallowWaterMultiLayerEquations{2,
                                                                                                         4,
                                                                                                         1},
                                                      dg, cache)
    (; boundary_condition_types, boundary_indices) = boundary_conditions
    (; contravariant_vectors) = cache.elements

    (; boundaries) = cache
    index_range = eachnode(dg)

    Trixi.foreach_enumerate(boundary_condition_types) do (i, boundary_condition)
        for boundary in boundary_indices[i]
            element = boundaries.neighbor_ids[boundary]
            node_indices = boundaries.node_indices[boundary]
            direction = Trixi.indices2direction(node_indices)

            i_node_start, i_node_step = Trixi.index_to_start_step_2d(node_indices[1],
                                                                     index_range)
            j_node_start, j_node_step = Trixi.index_to_start_step_2d(node_indices[2],
                                                                     index_range)

            i_node = i_node_start
            j_node = j_node_start
            for i in eachnode(dg)
                normal_direction = Trixi.get_normal_direction(direction,
                                                              contravariant_vectors,
                                                              i_node, j_node, element)

                u_inner = get_node_vars(u, equations, dg, i_node, j_node, element)

                u_outer = Trixi.get_boundary_outer_state(u_inner, t, boundary_condition,
                                                         normal_direction,
                                                         mesh, equations, dg, cache,
                                                         i_node, j_node, element)
                var_outer = u_outer[variable]

                # If the water height is specified as the limiting variable, apply limiting with the
                # total water height `H = h + b` instead.
                if variable == 1
                    var_outer = u_outer[variable] + u_outer[end]
                end

                var_min[i_node, j_node, element] = min(var_min[i_node, j_node, element],
                                                       var_outer)
                var_max[i_node, j_node, element] = max(var_max[i_node, j_node, element],
                                                       var_outer)

                i_node += i_node_step
                j_node += j_node_step
            end
        end
    end

    return nothing
end
