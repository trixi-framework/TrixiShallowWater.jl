# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# TODO: once working the mortar methods could likely be extended to the other
# equations types available in the package. Although for the multilayer equations
# care must be taken beacuse the pressure term is separated from the physical flux
# and directly placed in the nonconservative flux

# The methods below are specialized on the mortar type
# and called from the basic `create_cache` method at the top.
# TODO: better comment why we need more storage
function Trixi.create_cache(mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                            equations::ShallowWaterEquationsWetDry2D,
                            mortar_l2::Trixi.LobattoLegendreMortarL2, uEltype)
    # TODO: Taal performance using different types
    MA2d = Trixi.MArray{Tuple{nvariables(equations), nnodes(mortar_l2)},
                        uEltype, 2,
                        nvariables(equations) * nnodes(mortar_l2)}
    fstar_primary_upper_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    fstar_primary_lower_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    fstar_secondary_upper_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    fstar_secondary_lower_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    u_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]

    (; fstar_primary_upper_threaded, fstar_primary_lower_threaded,
       fstar_secondary_upper_threaded, fstar_secondary_lower_threaded, u_threaded)
end

function Trixi.prolong2mortars!(cache, u,
                                mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                equations::ShallowWaterEquationsWetDry2D,
                                mortar_l2::Trixi.LobattoLegendreMortarL2,
                                surface_integral, dg::DGSEM)
    @unpack neighbor_ids, node_indices = cache.mortars
    index_range = eachnode(dg)

    Trixi.@threaded for mortar in Trixi.eachmortar(dg, cache)
        # Copy solution data from the small elements using "delayed indexing" with
        # a start value and a step size to get the correct face and orientation.
        small_indices = node_indices[1, mortar]

        i_small_start, i_small_step = Trixi.index_to_start_step_2d(small_indices[1],
                                                                   index_range)
        j_small_start, j_small_step = Trixi.index_to_start_step_2d(small_indices[2],
                                                                   index_range)

        for position in 1:2
            i_small = i_small_start
            j_small = j_small_start
            element = neighbor_ids[position, mortar]
            for i in eachnode(dg)
                for v in eachvariable(equations)
                    cache.mortars.u[1, v, position, i, mortar] = u[v, i_small, j_small,
                                                                   element]
                end
                i_small += i_small_step
                j_small += j_small_step
            end
        end

        # Buffer to copy solution values of the large element in the correct orientation
        # before interpolating
        u_buffer = cache.u_threaded[Threads.threadid()]

        # Copy solution of large element face to buffer in the
        # correct orientation
        large_indices = node_indices[2, mortar]

        i_large_start, i_large_step = Trixi.index_to_start_step_2d(large_indices[1],
                                                                   index_range)
        j_large_start, j_large_step = Trixi.index_to_start_step_2d(large_indices[2],
                                                                   index_range)

        i_large = i_large_start
        j_large = j_large_start
        element = neighbor_ids[3, mortar]
        for i in eachnode(dg)
            # Compute and save the sigma variable from Benov et al. (essentially H = h+b),
            # momenta, and bottom topography into the buffer for projection
            temp = zeros(eltype(cache.mortars.u), 4)
            if u[1, i_large, j_large, element] > equations.threshold_wet
                temp[1] = u[1, i_large, j_large, element] + u[4, i_large, j_large, element]
            else
                temp[1] = equations.H0 # b + threshold?
            end
            temp[2:4] = u[2:4, i_large, j_large, element]

            for v in eachvariable(equations)
                u_buffer[v, i] = temp[v]
                # u_buffer[v, i] = u[v, i_large, j_large, element]
            end
            i_large += i_large_step
            j_large += j_large_step
        end

        # Interpolate large element face data from buffer to small face locations
        Trixi.multiply_dimensionwise!(view(cache.mortars.u, 2, :, 1, :, mortar),
                                      mortar_l2.forward_lower,
                                      u_buffer)
        Trixi.multiply_dimensionwise!(view(cache.mortars.u, 2, :, 2, :, mortar),
                                      mortar_l2.forward_upper,
                                      u_buffer)
    end

    return nothing
end

function Trixi.calc_mortar_flux!(surface_flux_values,
                                 mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                 nonconservative_terms,
                                 equations::ShallowWaterEquationsWetDry2D,
                                 mortar_l2::Trixi.LobattoLegendreMortarL2,
                                 surface_integral, dg::DG, cache)
    @unpack neighbor_ids, node_indices = cache.mortars
    @unpack contravariant_vectors = cache.elements
    @unpack fstar_primary_upper_threaded, fstar_primary_lower_threaded, fstar_secondary_upper_threaded, fstar_secondary_lower_threaded = cache
    surface_flux, nonconservative_flux = surface_integral.surface_flux
    index_range = eachnode(dg)

    Trixi.@threaded for mortar in Trixi.eachmortar(dg, cache)
        # Choose thread-specific pre-allocated containers
        # for the numerical flux on the small elements as well as
        # the physical flux evaluated at the projected solution
        # from the large element
        fstar_primary = (fstar_primary_lower_threaded[Threads.threadid()],
                         fstar_primary_upper_threaded[Threads.threadid()])

        fstar_secondary = (fstar_secondary_lower_threaded[Threads.threadid()],
                           fstar_secondary_upper_threaded[Threads.threadid()])

        # Get index information on the small elements
        small_indices = node_indices[1, mortar]
        small_direction = Trixi.indices2direction(small_indices)

        i_small_start, i_small_step = Trixi.index_to_start_step_2d(small_indices[1],
                                                                   index_range)
        j_small_start, j_small_step = Trixi.index_to_start_step_2d(small_indices[2],
                                                                   index_range)

        for position in 1:2
            i_small = i_small_start
            j_small = j_small_start
            element = neighbor_ids[position, mortar]
            for node in eachnode(dg)
                # Get the normal direction on the small element.
                # Note, contravariant vectors at interfaces in negative coordinate direction
                # are pointing inwards. This is handled by `get_normal_direction`.
                normal_direction = Trixi.get_normal_direction(small_direction,
                                                              contravariant_vectors,
                                                              i_small, j_small, element)

                # Compute the numerical flux in the normal direction on the small elements
                Trixi.calc_mortar_flux!(fstar_primary, fstar_secondary,
                                        mesh, nonconservative_terms, equations,
                                        surface_integral, dg, cache,
                                        mortar, position, normal_direction,
                                        node)

                i_small += i_small_step
                j_small += j_small_step
            end
        end

        # Buffer to interpolate flux values of the large element to before
        # copying in the correct orientation
        u_buffer = cache.u_threaded[Threads.threadid()]

        # in calc_interface_flux!, the interface flux is computed once over each
        # interface using the normal from the "primary" element. The result is then
        # passed back to the "secondary" element, flipping the sign to account for the
        # change in the normal direction. For mortars, this sign flip occurs in
        # "mortar_fluxes_to_elements!" instead.
        Trixi.mortar_fluxes_to_elements!(surface_flux_values,
                                         mesh, equations, mortar_l2, dg, cache,
                                         mortar, fstar_primary, fstar_secondary, u_buffer)
    end

    return nothing
end

# Inlined version of the mortar flux computation on small elements for equations with conservative and
# nonconservative terms
@inline function Trixi.calc_mortar_flux!(fstar_primary, fstar_secondary,
                                         mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                         nonconservative_terms::True,
                                         equations::ShallowWaterEquationsWetDry2D,
                                         surface_integral, dg::DG, cache,
                                         mortar_index, position_index, normal_direction,
                                         node_index)
    @unpack u = cache.mortars
    surface_flux, nonconservative_flux = surface_integral.surface_flux

    u_ll, u_rr = Trixi.get_surface_node_vars(u, equations, dg, position_index, node_index,
                                             mortar_index)
    # TODO: The surface values in `u_rr` are still in terms of sigma, momenta, and bottom topography
    # (1) unpack the sigma to create the water height on the mortar from Eq. 41 in Benov et al.
    # (2) perform hydrostatic reconstruction on the mortar solution
    # (3) use these HR quantities to compute the flux and nonconservative terms
    # These last two steps, as far as I recall, occur within the surface flux from the
    # provided solution states
    u_rr_temp = zeros(eltype(u_rr), 4)

    u_rr_temp[1] = max(u_rr[1] - u_rr[4], equations.threshold_wet) # zero(eltype(u_rr)))
    u_rr_temp[2:4] = u_rr[2:4]

    # Compute conservative flux
    flux = surface_flux(u_ll, u_rr_temp, normal_direction, equations)
    # flux = surface_flux(u_ll, u_rr, normal_direction, equations)

    # Compute nonconservative flux and add it to the conservative flux.
    # The nonconservative flux is scaled by a factor of 0.5 based on
    # the interpretation of global SBP operators coupled discontinuously via
    # central fluxes/SATs
    noncons_primary = nonconservative_flux(u_ll, u_rr_temp,
                                           normal_direction, normal_direction,
                                           equations)
    noncons_secondary = nonconservative_flux(u_rr_temp, u_ll,
                                             normal_direction, normal_direction,
                                             equations)
    # noncons_primary = nonconservative_flux(u_ll, u_rr,
    #                                        normal_direction, normal_direction,
    #                                        equations)
    # noncons_secondary = nonconservative_flux(u_rr, u_ll,
    #                                          normal_direction, normal_direction,
    #                                          equations)

    flux_plus_noncons_primary = flux + 0.5 * noncons_primary
    flux_plus_noncons_secondary = flux + 0.5 * noncons_secondary

    # Copy to buffer
    set_node_vars!(fstar_primary[position_index], flux_plus_noncons_primary, equations, dg, node_index)
    set_node_vars!(fstar_secondary[position_index], flux_plus_noncons_secondary, equations, dg, node_index)
end

@inline function Trixi.mortar_fluxes_to_elements!(surface_flux_values,
                                                  mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                                  equations::ShallowWaterEquationsWetDry2D,
                                                  mortar_l2::Trixi.LobattoLegendreMortarL2,
                                                  dg::DGSEM, cache, mortar,
                                                  fstar_primary, fstar_secondary, u_buffer)
    @unpack contravariant_vectors = cache.elements
    @unpack neighbor_ids, node_indices = cache.mortars
    surface_flux, nonconservative_flux = dg.surface_integral.surface_flux

    # Copy solution small to small
    small_indices = node_indices[1, mortar]
    small_direction = Trixi.indices2direction(small_indices)

    for position in 1:2
        element = neighbor_ids[position, mortar]
        for i in eachnode(dg)
            for v in eachvariable(equations)
                surface_flux_values[v, i, small_direction, element] = fstar_primary[position][v,
                                                                                      i]
            end
        end
    end

    # Project small numerical fluxes and physical flux computed on the projected
    # large element solution back onto large element.
    # This is bascially Eq. (46) from Benov et al. where the factor of
    # 1/2 is already included in `reverse_upper` and `reverse_lower`.
    Trixi.multiply_dimensionwise!(u_buffer,
                            # mortar_l2.reverse_upper, fstar_primary[2],
                            # mortar_l2.reverse_lower, fstar_primary[1])
                            mortar_l2.reverse_upper, fstar_secondary[2],
                            mortar_l2.reverse_lower, fstar_secondary[1])

    # The flux is calculated in the outward direction of the small elements,
    # so the sign must be switched to get the flux in outward direction
    # of the large element.
    # The contravariant vectors of the large element (and therefore the normal
    # vectors of the large element as well) are twice as large as the
    # contravariant vectors of the small elements. Therefore, the flux needs
    # to be scaled by a factor of 2 to obtain the flux of the large element.
    u_buffer .*= -2

    # Copy interpolated flux values from buffer to large element face in the
    # correct orientation.
    # Note that the index of the small sides will always run forward but
    # the index of the large side might need to run backwards for flipped sides.
    large_element = neighbor_ids[3, mortar]
    large_indices = node_indices[2, mortar]
    large_direction = Trixi.indices2direction(large_indices)

    if :i_backward in large_indices
        for i in eachnode(dg)
            for v in eachvariable(equations)
                surface_flux_values[v, end + 1 - i, large_direction, large_element] = u_buffer[v, i]
            end
        end
    else
        for i in eachnode(dg)
            for v in eachvariable(equations)
                surface_flux_values[v, i, large_direction, large_element] = u_buffer[v, i]
            end
        end
    end

    return nothing
end
end # @muladd
