# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# TODO: Once working the mortar methods could likely be extended to the other
# equation types available in the package. Although for the multilayer equations
# care must be taken because the pressure term is separated from the physical flux
# and directly placed in the nonconservative flux

# The methods below are specialized on the `P4estShallowWaterMortarContainer`
# mortar type. Extra storage is necessary for the numerical flux plus nonconservative term
# on either side of the parent (large) element mortars. It also requires a work
# array to compute the flux with the projected large element solutions on each mortar
# to ensure we project back the flux penalty.
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
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
    f_upper_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    f_lower_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    u_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]
    f_threaded = MA2d[MA2d(undef) for _ in 1:Threads.nthreads()]

    (; fstar_primary_upper_threaded, fstar_primary_lower_threaded,
     fstar_secondary_upper_threaded, fstar_secondary_lower_threaded,
     f_upper_threaded, f_lower_threaded, u_threaded, f_threaded)
end

function Trixi.prolong2mortars!(cache, u,
                                mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                equations::ShallowWaterEquationsWetDry2D,
                                mortar_l2::Trixi.LobattoLegendreMortarL2,
                                dg::DGSEM)
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
        # before projection
        u_buffer = cache.u_threaded[Threads.threadid()]

        # Copy solution of large element face to buffer in the correct orientation
        large_indices = node_indices[2, mortar]

        i_large_start, i_large_step = Trixi.index_to_start_step_2d(large_indices[1],
                                                                   index_range)
        j_large_start, j_large_step = Trixi.index_to_start_step_2d(large_indices[2],
                                                                   index_range)

        i_large = i_large_start
        j_large = j_large_start
        element = neighbor_ids[3, mortar]
        for i in eachnode(dg)
            # This strategy from Benov et al. (https://doi.org/10.1016/j.jcp.2018.02.008) assumes that
            # we know a constant background water height `H0` which we perturb around. This may be restrictive
            # in practice but a good place to start with the development. We may need to consider more
            # sophisticated positivity preserving projections of the solution like those found
            # in the ALE-DG community for the Euler equations with gravity to remove this assumption.
            # That is, we might be able to directly project the water height `h` instead while maintaining
            # important steady-state solution behavior.
            # Note, a small shift is required to ensure we catch water heights close to the threshold
            if u[1, i_large, j_large, element] >=
               2 * (equations.threshold_limiter + eps())
                u_buffer[1, i] = u[1, i_large, j_large, element] +
                                 u[4, i_large, j_large, element]
            else
                u_buffer[1, i] = equations.H0
            end
            for v in 2:4
                u_buffer[v, i] = u[v, i_large, j_large, element]
            end

            for v in eachvariable(equations)
                # Save a copy of the (unprojected) parent solution on the mortar
                cache.mortars.u_parent[v, i, mortar] = u[v, i_large, j_large, element]
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

        # After the projection of the constant solution we modify the values
        # in the first solution variable to no longer be the total water height H = h+b
        # and instead recover the conservative water height variable `h`.
        # Basically, unpacking the total water height variable to create the projected local water
        # height `h` from Eq. 41 in Benov et al.
        for i in eachnode(dg)
            cache.mortars.u[2, 1, 1, i, mortar] = max(cache.mortars.u[2, 1, 1, i,
                                                                      mortar] -
                                                      cache.mortars.u[2, 4, 1, i,
                                                                      mortar],
                                                      equations.threshold_limiter)
            cache.mortars.u[2, 1, 2, i, mortar] = max(cache.mortars.u[2, 1, 2, i,
                                                                      mortar] -
                                                      cache.mortars.u[2, 4, 2, i,
                                                                      mortar],
                                                      equations.threshold_limiter)

            # Safety application of velocity desingularization and water height cutoff on the mortars.
            # Here we use the rather conservative velocity desingularization tolerance of `1e-4`
            # which works for the AMR testing. The more standard value of `1e-6` works
            # for well-balancedness testing, but crashes the AMR simulation elixir.
            #
            # For details on the motivation of the velocity desingularization see
            # - A. Chertock, S. Cui, A. Kurganov, T. Wu (2015)
            #   Well-balanced positivity preserving central-upwind scheme for
            #   the shallow water system with friction terms
            #   [DOI: 10.1002/fld.4023](https://doi.org/10.1002/fld.4023)

            # Mortars with the copied or projected solution
            tol = 1e-4
            for child in 1:2, side in 1:2
                h = cache.mortars.u[side, 1, child, i, mortar]
                cache.mortars.u[side, 2, child, i, mortar] = h * (2 * h *
                                                              cache.mortars.u[side, 2,
                                                                              child, i,
                                                                              mortar]) /
                                                             (h^2 + max(h^2, tol))
                cache.mortars.u[side, 3, child, i, mortar] = h * (2 * h *
                                                              cache.mortars.u[side, 3,
                                                                              child, i,
                                                                              mortar]) /
                                                             (h^2 + max(h^2, tol))
                if cache.mortars.u[side, 1, child, i, mortar] <=
                   equations.threshold_limiter
                    cache.mortars.u[side, 1, child, i, mortar] = equations.threshold_limiter
                    cache.mortars.u[side, 2, child, i, mortar] = zero(eltype(u))
                    cache.mortars.u[side, 3, child, i, mortar] = zero(eltype(u))
                end
            end

            # Unprojected parent solution copied into the mortar storage
            h = cache.mortars.u_parent[1, i, mortar]
            cache.mortars.u_parent[2, i, mortar] = h * (2 * h *
                                                    cache.mortars.u_parent[2, i, mortar]) /
                                                   (h^2 + max(h^2, tol))
            cache.mortars.u_parent[3, i, mortar] = h * (2 * h *
                                                    cache.mortars.u_parent[3, i, mortar]) /
                                                   (h^2 + max(h^2, tol))
            if cache.mortars.u_parent[1, i, mortar] <= equations.threshold_limiter
                cache.mortars.u_parent[1, i, mortar] = equations.threshold_limiter
                cache.mortars.u_parent[2, i, mortar] = zero(eltype(u))
                cache.mortars.u_parent[3, i, mortar] = zero(eltype(u))
            end
        end
    end

    return nothing
end

# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function Trixi.calc_mortar_flux!(surface_flux_values,
                                 mesh::Union{P4estMesh{2}, T8codeMesh{2}},
                                 nonconservative_terms,
                                 equations::ShallowWaterEquationsWetDry2D,
                                 mortar_l2::Trixi.LobattoLegendreMortarL2,
                                 surface_integral, dg::DG, cache)
    @unpack neighbor_ids, node_indices = cache.mortars
    @unpack contravariant_vectors = cache.elements
    @unpack (fstar_primary_upper_threaded, fstar_primary_lower_threaded,
    fstar_secondary_upper_threaded, fstar_secondary_lower_threaded,
    f_upper_threaded, f_lower_threaded) = cache

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

        f = (f_lower_threaded[Threads.threadid()],
             f_upper_threaded[Threads.threadid()])

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

                # The projected solution from the large element is always stored in `u_rr`
                _, u_rr = Trixi.get_surface_node_vars(cache.mortars.u, equations, dg,
                                                      position, node, mortar)

                # Compute conservative flux. Note, we use the surface flux here evaluated at the same
                # solution state to recover the physical flux at this point because the surface flux
                # has in-built mechanisms to avoid division by zero in dry regions whereas `Trixi.flux`
                # does not have such mechanisms to desingularize the velocity computation.
                flux = surface_flux(u_rr, u_rr, normal_direction, equations)

                # Compute nonconservative flux and add it to the conservative flux.
                # The nonconservative flux is scaled by a factor of 0.5 based on
                # the interpretation of global SBP operators coupled discontinuously via
                # central fluxes/SATs
                noncons = nonconservative_flux(u_rr, u_rr, normal_direction,
                                               equations)

                flux_plus_noncons = flux + 0.5f0 * noncons

                # Copy to the physical flux buffer
                set_node_vars!(f[position], flux_plus_noncons, equations, dg, node)

                i_small += i_small_step
                j_small += j_small_step
            end
        end

        # Buffer to interpolate flux values of the large element to before
        # copying in the correct orientation
        u_buffer = cache.u_threaded[Threads.threadid()]

        # In calc_interface_flux!, the interface flux is computed once over each
        # interface using the normal from the "primary" element. The result is then
        # passed back to the "secondary" element, flipping the sign to account for the
        # change in the normal direction. For mortars, this sign flip occurs in
        # "mortar_fluxes_to_elements!" instead.
        Trixi.mortar_fluxes_to_elements!(surface_flux_values,
                                         mesh, equations, mortar_l2, dg, cache,
                                         mortar, fstar_primary, fstar_secondary,
                                         f, u_buffer)
    end

    return nothing
end

# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
@inline function Trixi.mortar_fluxes_to_elements!(surface_flux_values,
                                                  mesh::Union{P4estMesh{2},
                                                              T8codeMesh{2}},
                                                  equations::ShallowWaterEquationsWetDry2D,
                                                  mortar_l2::Trixi.LobattoLegendreMortarL2,
                                                  dg::DGSEM, cache, mortar,
                                                  fstar_primary, fstar_secondary,
                                                  f_large, u_buffer)
    @unpack contravariant_vectors = cache.elements
    @unpack neighbor_ids, node_indices = cache.mortars
    surface_flux, nonconservative_flux = dg.surface_integral.surface_flux

    # Copy surface flux data from small to small
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

    # Project small numerical fluxes and physical flux penalty computed on the projected
    # large element solution back onto large element.
    # This is basically Eq. (46) from Benov et al. where the factor of 1/2 is already
    # already included in `reverse_upper` and `reverse_lower` operators.
    penalty_upper = Trixi.SMatrix(fstar_secondary[2] - f_large[2])
    penalty_lower = Trixi.SMatrix(fstar_secondary[1] - f_large[1])
    Trixi.multiply_dimensionwise!(u_buffer,
                                  mortar_l2.reverse_upper, penalty_upper,
                                  mortar_l2.reverse_lower, penalty_lower)

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

    # From the unprojected solution stored in the mortars we have access to compute the flux
    # on the parent elements and remove the physical flux evaluated
    # at the unprojected solution state that is present from the volume integral computation
    # later in the computation.
    index_range = eachnode(dg)
    i_large_start, i_large_step = Trixi.index_to_start_step_2d(large_indices[1],
                                                               index_range)
    j_large_start, j_large_step = Trixi.index_to_start_step_2d(large_indices[2],
                                                               index_range)

    i_large = i_large_start
    j_large = j_large_start
    flux_buffer = cache.f_threaded[Threads.threadid()]
    for node in eachnode(dg)
        # Get the proper normal_direction now that we are back computing on the large element
        normal_direction = Trixi.get_normal_direction(large_direction,
                                                      contravariant_vectors,
                                                      i_large, j_large, large_element)

        # Compute conservative flux. Note, we use the surface flux here evaluated at the same
        # solution state to recover the physical flux at this point because the surface flux
        # has in-built mechanisms to avoid division by zero in dry regions whereas `Trixi.flux`
        # does not have such mechanisms to desingularize the velocity computation.
        flux = surface_flux(view(cache.mortars.u_parent, :, node, mortar),
                            view(cache.mortars.u_parent, :, node, mortar),
                            normal_direction, equations)

        noncons = nonconservative_flux(view(cache.mortars.u_parent, :, node, mortar),
                                       view(cache.mortars.u_parent, :, node, mortar),
                                       normal_direction,
                                       equations)

        flux_plus_noncons = flux + 0.5f0 * noncons

        set_node_vars!(flux_buffer, flux_plus_noncons, equations, dg, node)

        i_large += i_large_step
        j_large += j_large_step
    end

    # Compute the surface flux values on the large element. Note, we have projected
    # the flux penalty from the small element mortars back onto the large element side.
    # However, the large element has a local contribution to the strong form flux penalty
    # that must be removed to ensure consistency. This is done via the values precomputed
    # and stored above in the `flux_buffer` which is the physical flux and nonconservative
    # term evaluated at the unprojected large element solution and `normal_direction`.
    if :i_backward in large_indices
        for i in eachnode(dg)
            for v in eachvariable(equations)
                surface_flux_values[v, end + 1 - i, large_direction, large_element] = u_buffer[v,
                                                                                               i] +
                                                                                      flux_buffer[v,
                                                                                                  i]
            end
        end
    else
        for i in eachnode(dg)
            for v in eachvariable(equations)
                surface_flux_values[v, i, large_direction, large_element] = u_buffer[v,
                                                                                     i] +
                                                                            flux_buffer[v,
                                                                                        i]
            end
        end
    end

    return nothing
end
end # @muladd
