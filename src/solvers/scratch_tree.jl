# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

function prolong2mortars!(cache, u,
                          mesh::TreeMesh{2}, equations,
                          mortar_l2::LobattoLegendreMortarL2, surface_integral,
                          dg::DGSEM)
    @threaded for mortar in eachmortar(dg, cache)
        large_element = cache.mortars.neighbor_ids[3, mortar]
        upper_element = cache.mortars.neighbor_ids[2, mortar]
        lower_element = cache.mortars.neighbor_ids[1, mortar]

        # Copy solution small to small
        if cache.mortars.large_sides[mortar] == 1 # -> small elements on right side
            if cache.mortars.orientations[mortar] == 1
                # L2 mortars in x-direction
                for l in eachnode(dg)
                    for v in eachvariable(equations)
                        cache.mortars.u_upper[2, v, l, mortar] = u[v, 1, l,
                                                                   upper_element]
                        cache.mortars.u_lower[2, v, l, mortar] = u[v, 1, l,
                                                                   lower_element]
                    end
                end
            else
                # L2 mortars in y-direction
                for l in eachnode(dg)
                    for v in eachvariable(equations)
                        cache.mortars.u_upper[2, v, l, mortar] = u[v, l, 1,
                                                                   upper_element]
                        cache.mortars.u_lower[2, v, l, mortar] = u[v, l, 1,
                                                                   lower_element]
                    end
                end
            end
        else # large_sides[mortar] == 2 -> small elements on left side
            if cache.mortars.orientations[mortar] == 1
                # L2 mortars in x-direction
                for l in eachnode(dg)
                    for v in eachvariable(equations)
                        cache.mortars.u_upper[1, v, l, mortar] = u[v, nnodes(dg), l,
                                                                   upper_element]
                        cache.mortars.u_lower[1, v, l, mortar] = u[v, nnodes(dg), l,
                                                                   lower_element]
                    end
                end
            else
                # L2 mortars in y-direction
                for l in eachnode(dg)
                    for v in eachvariable(equations)
                        cache.mortars.u_upper[1, v, l, mortar] = u[v, l, nnodes(dg),
                                                                   upper_element]
                        cache.mortars.u_lower[1, v, l, mortar] = u[v, l, nnodes(dg),
                                                                   lower_element]
                    end
                end
            end
        end

        # Interpolate large element face data to small interface locations
        if cache.mortars.large_sides[mortar] == 1 # -> large element on left side
            leftright = 1
            if cache.mortars.orientations[mortar] == 1
                # L2 mortars in x-direction
                u_large = view(u, :, nnodes(dg), :, large_element)
                element_solutions_to_mortars!(cache.mortars, mortar_l2, leftright,
                                              mortar, u_large)
            else
                # L2 mortars in y-direction
                u_large = view(u, :, :, nnodes(dg), large_element)
                element_solutions_to_mortars!(cache.mortars, mortar_l2, leftright,
                                              mortar, u_large)
            end
        else # large_sides[mortar] == 2 -> large element on right side
            leftright = 2
            if cache.mortars.orientations[mortar] == 1
                # L2 mortars in x-direction
                u_large = view(u, :, 1, :, large_element)
                element_solutions_to_mortars!(cache.mortars, mortar_l2, leftright,
                                              mortar, u_large)
            else
                # L2 mortars in y-direction
                u_large = view(u, :, :, 1, large_element)
                element_solutions_to_mortars!(cache.mortars, mortar_l2, leftright,
                                              mortar, u_large)
            end
        end
    end

    return nothing
end

@inline function element_solutions_to_mortars!(mortars,
                                               mortar_l2::LobattoLegendreMortarL2,
                                               leftright, mortar,
                                               u_large::AbstractArray{<:Any, 2})
    multiply_dimensionwise!(view(mortars.u_upper, leftright, :, :, mortar),
                            mortar_l2.forward_upper, u_large)
    multiply_dimensionwise!(view(mortars.u_lower, leftright, :, :, mortar),
                            mortar_l2.forward_lower, u_large)
    return nothing
end

function calc_mortar_flux!(surface_flux_values,
                           mesh::TreeMesh{2},
                           nonconservative_terms::True, equations,
                           mortar_l2::LobattoLegendreMortarL2,
                           surface_integral, dg::DG, cache)
    surface_flux, nonconservative_flux = surface_integral.surface_flux
    @unpack u_lower, u_upper, orientations, large_sides = cache.mortars
    @unpack fstar_upper_threaded, fstar_lower_threaded = cache

    @threaded for mortar in eachmortar(dg, cache)
        # Choose thread-specific pre-allocated container
        fstar_upper = fstar_upper_threaded[Threads.threadid()]
        fstar_lower = fstar_lower_threaded[Threads.threadid()]

        # Calculate fluxes
        orientation = orientations[mortar]
        calc_fstar!(fstar_upper, equations, surface_flux, dg, u_upper, mortar,
                    orientation)
        calc_fstar!(fstar_lower, equations, surface_flux, dg, u_lower, mortar,
                    orientation)

        # Add nonconservative fluxes.
        # These need to be adapted on the geometry (left/right) since the order of
        # the arguments matters, based on the global SBP operator interpretation.
        # The same interpretation (global SBP operators coupled discontinuously via
        # central fluxes/SATs) explains why we need the factor 0.5.
        # Alternatively, you can also follow the argumentation of Bohm et al. 2018
        # ("nonconservative diamond flux")
        if large_sides[mortar] == 1 # -> small elements on right side
            for i in eachnode(dg)
                # Pull the left and right solutions
                u_upper_ll, u_upper_rr = get_surface_node_vars(u_upper, equations, dg,
                                                               i, mortar)
                u_lower_ll, u_lower_rr = get_surface_node_vars(u_lower, equations, dg,
                                                               i, mortar)
                # Call pointwise nonconservative term
                noncons_upper = nonconservative_flux(u_upper_ll, u_upper_rr,
                                                     orientation, equations)
                noncons_lower = nonconservative_flux(u_lower_ll, u_lower_rr,
                                                     orientation, equations)
                # Add to primary and secondary temporary storage
                multiply_add_to_node_vars!(fstar_upper, 0.5, noncons_upper, equations,
                                           dg, i)
                multiply_add_to_node_vars!(fstar_lower, 0.5, noncons_lower, equations,
                                           dg, i)
            end
        else # large_sides[mortar] == 2 -> small elements on the left
            for i in eachnode(dg)
                # Pull the left and right solutions
                u_upper_ll, u_upper_rr = get_surface_node_vars(u_upper, equations, dg,
                                                               i, mortar)
                u_lower_ll, u_lower_rr = get_surface_node_vars(u_lower, equations, dg,
                                                               i, mortar)
                # Call pointwise nonconservative term
                noncons_upper = nonconservative_flux(u_upper_rr, u_upper_ll,
                                                     orientation, equations)
                noncons_lower = nonconservative_flux(u_lower_rr, u_lower_ll,
                                                     orientation, equations)
                # Add to primary and secondary temporary storage
                multiply_add_to_node_vars!(fstar_upper, 0.5, noncons_upper, equations,
                                           dg, i)
                multiply_add_to_node_vars!(fstar_lower, 0.5, noncons_lower, equations,
                                           dg, i)
            end
        end

        mortar_fluxes_to_elements!(surface_flux_values,
                                   mesh, equations, mortar_l2, dg, cache,
                                   mortar, fstar_upper, fstar_lower)
    end

    return nothing
end

@inline function calc_fstar!(destination::AbstractArray{<:Any, 2}, equations,
                             surface_flux, dg::DGSEM,
                             u_interfaces, interface, orientation)
    for i in eachnode(dg)
        # Call pointwise two-point numerical flux function
        u_ll, u_rr = get_surface_node_vars(u_interfaces, equations, dg, i, interface)
        flux = surface_flux(u_ll, u_rr, orientation, equations)

        # Copy flux to left and right element storage
        set_node_vars!(destination, flux, equations, dg, i)
    end

    return nothing
end

@inline function mortar_fluxes_to_elements!(surface_flux_values,
                                            mesh::TreeMesh{2}, equations,
                                            mortar_l2::LobattoLegendreMortarL2,
                                            dg::DGSEM, cache,
                                            mortar, fstar_upper, fstar_lower)
    large_element = cache.mortars.neighbor_ids[3, mortar]
    upper_element = cache.mortars.neighbor_ids[2, mortar]
    lower_element = cache.mortars.neighbor_ids[1, mortar]

    # Copy flux small to small
    if cache.mortars.large_sides[mortar] == 1 # -> small elements on right side
        if cache.mortars.orientations[mortar] == 1
            # L2 mortars in x-direction
            direction = 1
        else
            # L2 mortars in y-direction
            direction = 3
        end
    else # large_sides[mortar] == 2 -> small elements on left side
        if cache.mortars.orientations[mortar] == 1
            # L2 mortars in x-direction
            direction = 2
        else
            # L2 mortars in y-direction
            direction = 4
        end
    end
    surface_flux_values[:, :, direction, upper_element] .= fstar_upper
    surface_flux_values[:, :, direction, lower_element] .= fstar_lower

    # Project small fluxes to large element
    if cache.mortars.large_sides[mortar] == 1 # -> large element on left side
        if cache.mortars.orientations[mortar] == 1
            # L2 mortars in x-direction
            direction = 2
        else
            # L2 mortars in y-direction
            direction = 4
        end
    else # large_sides[mortar] == 2 -> large element on right side
        if cache.mortars.orientations[mortar] == 1
            # L2 mortars in x-direction
            direction = 1
        else
            # L2 mortars in y-direction
            direction = 3
        end
    end

    # TODO: Taal performance
    # for v in eachvariable(equations)
    #   # The code below is semantically equivalent to
    #   # surface_flux_values[v, :, direction, large_element] .=
    #   #   (mortar_l2.reverse_upper * fstar_upper[v, :] + mortar_l2.reverse_lower * fstar_lower[v, :])
    #   # but faster and does not allocate.
    #   # Note that `true * some_float == some_float` in Julia, i.e. `true` acts as
    #   # a universal `one`. Hence, the second `mul!` means "add the matrix-vector
    #   # product to the current value of the destination".
    #   @views mul!(surface_flux_values[v, :, direction, large_element],
    #               mortar_l2.reverse_upper, fstar_upper[v, :])
    #   @views mul!(surface_flux_values[v, :, direction, large_element],
    #               mortar_l2.reverse_lower,  fstar_lower[v, :], true, true)
    # end
    # The code above could be replaced by the following code. However, the relative efficiency
    # depends on the types of fstar_upper/fstar_lower and dg.l2mortar_reverse_upper.
    # Using StaticArrays for both makes the code above faster for common test cases.
    multiply_dimensionwise!(view(surface_flux_values, :, :, direction, large_element),
                            mortar_l2.reverse_upper, fstar_upper,
                            mortar_l2.reverse_lower, fstar_lower)

    return nothing
end
end # @muladd
