# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Refine elements in the DG solver based on a list of cell_ids that should be refined.
# On `P4estMesh`, if an element refines the solution scaled by the Jacobian `J*u` is interpolated
# from the parent element into the four children elements. The solution on each child
# element is then recovered by dividing by the new element Jacobians.
# After refinement, projections may introduce inadmissible solution states like negative
# water heights. So the limiter is called to remedy this.
#
# TODO: Requires modification for use with `TreeMesh` because the Jacobian is constant.
# A generic version is available in Trixi.jl in `src/callbacks_step/amr_dg2d.jl`.
# However, the well-balanced mortar implementation would need added to `TreeMesh` first.
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function Trixi.refine!(u_ode::AbstractVector, adaptor, mesh::P4estMesh{2},
                       equations::ShallowWaterEquationsWetDry2D,
                       dg::DGSEM, cache, elements_to_refine)
    # Return early if there is nothing to do
    if isempty(elements_to_refine)
        # TODO: MPI parallel runs unavailable. Requires analogous implementations
        # of `calc_mpi_interface_flux!` and `calc_mpi_mortar_flux!` in Trixi.jl
        # `src/solvers/dgsem_p4est/dg_2d_parallel.jl`.
        # if Trixi.mpi_isparallel()
        #     # MPICache init uses all-to-all communication -> reinitialize even if there is nothing to do
        #     # locally (there still might be other MPI ranks that have refined elements)
        #     Trixi.reinitialize_containers!(mesh, equations, dg, cache)
        # end
        return
    end

    # Determine for each existing element whether it needs to be refined
    needs_refinement = falses(nelements(dg, cache))
    needs_refinement[elements_to_refine] .= true

    # Retain current solution data
    old_n_elements = nelements(dg, cache)
    old_u_ode = copy(u_ode)
    old_inverse_jacobian = copy(cache.elements.inverse_jacobian)
    # OBS! If we don't GC.@preserve old_u_ode and old_inverse_jacobian, they might be GC'ed
    GC.@preserve old_u_ode old_inverse_jacobian begin
        old_u = Trixi.wrap_array(old_u_ode, mesh, equations, dg, cache)

        # Loop over all elements in old container and scale the old solution by the Jacobian
        # prior to projection
        for old_element_id in 1:old_n_elements
            for v in eachvariable(equations)
                old_u[v, .., old_element_id] .= (old_u[v, .., old_element_id] ./
                                                 old_inverse_jacobian[..,
                                                                      old_element_id])
            end
        end

        Trixi.reinitialize_containers!(mesh, equations, dg, cache)

        resize!(u_ode,
                nvariables(equations) * nnodes(dg)^ndims(mesh) * nelements(dg, cache))
        u = Trixi.wrap_array(u_ode, mesh, equations, dg, cache)

        # Loop over all elements in old container and either copy them or refine them
        element_id = 1
        for old_element_id in 1:old_n_elements
            if needs_refinement[old_element_id]
                # Refine element and store solution directly in new data structure
                Trixi.refine_element!(u, element_id, old_u, old_element_id,
                                      adaptor, equations, dg)

                # Before `element_id` is incremented, divide by the new Jacobians on each
                # child element and save the result
                for m in 0:3 # loop over the children
                    for v in eachvariable(equations)
                        u[v, .., element_id + m] .*= (0.25f0 .*
                                                      cache.elements.inverse_jacobian[..,
                                                                                      element_id + m])
                    end
                end

                # Increment `element_id` on the refined mesh with the number
                # of children, i.e., 4 in 2D
                element_id += 2^ndims(mesh)
            else
                # Copy old element data to new element container and remove Jacobian scaling
                for v in eachvariable(equations)
                    u[v, .., element_id] .= (old_u[v, .., old_element_id] .*
                                             old_inverse_jacobian[..,
                                                                  old_element_id])
                end

                # No refinement occurred, so increment `element_id` on the new mesh by one
                element_id += 1
            end
        end
        # If everything is correct, we should have processed all elements.
        # Depending on whether the last element processed above had to be refined or not,
        # the counter `element_id` can have two different values at the end.
        @assert element_id ==
                nelements(dg, cache) +
                1||element_id == nelements(dg, cache) + 2^ndims(mesh) "element_id = $element_id, nelements(dg, cache) = $(nelements(dg, cache))"
    end # GC.@preserve old_u_ode old_inverse_jacobian

    # Apply limiter to ensure admissible solution states
    limiter_shallow_water!(u, equations.threshold_limiter, waterheight,
                           mesh, equations, dg, cache)
    return nothing
end

# Coarsen elements in the DG solver based on a list of cell_ids that should be removed.
# On `P4estMesh`, if an element coarsens the solution scaled by the Jacobian `J*u` is projected
# from the four children elements back onto the parent element. The solution on the parent
# element is then recovered by dividing by the new element Jacobian.
# After coarsening, projections may introduce inadmissible solution states like negative
# water heights. So the limiter is called to remedy this.
#
# TODO: Requires modification for use with `TreeMesh` because the Jacobian is constant.
# A generic version is available in Trixi.jl in `src/callbacks_step/amr_dg2d.jl`
# However, the well-balanced mortar implementation would need added to `TreeMesh` first.
#
# !!! warning "Experimental code"
#     This is an experimental feature and may change in future releases.
function Trixi.coarsen!(u_ode::AbstractVector, adaptor, mesh::P4estMesh{2},
                        equations::ShallowWaterEquationsWetDry2D,
                        dg::DGSEM, cache, elements_to_remove)
    # Return early if there is nothing to do
    if isempty(elements_to_remove)
        # TODO: MPI parallel runs unavailable. Requires analogous implementations
        # of `calc_mpi_interface_flux!` and `calc_mpi_mortar_flux!` in Trixi.jl
        # `src/solvers/dgsem_p4est/dg_2d_parallel.jl`.
        # if Trixi.mpi_isparallel()
        #     # MPICache init uses all-to-all communication -> reinitialize even if there is nothing to do
        #     # locally (there still might be other MPI ranks that have coarsened elements)
        #     Trixi.reinitialize_containers!(mesh, equations, dg, cache)
        # end
        return
    end

    # Determine for each old element whether it needs to be removed
    to_be_removed = falses(nelements(dg, cache))
    to_be_removed[elements_to_remove] .= true

    # Retain current solution data and Jacobians
    old_n_elements = nelements(dg, cache)
    old_u_ode = copy(u_ode)
    old_inverse_jacobian = copy(cache.elements.inverse_jacobian)
    # OBS! If we don't GC.@preserve old_u_ode and old_inverse_jacobian, they might be GC'ed
    GC.@preserve old_u_ode old_inverse_jacobian begin
        old_u = Trixi.wrap_array(old_u_ode, mesh, equations, dg, cache)

        # Loop over all elements in old container and scale the old solution by the Jacobian
        # prior to projection
        for old_element_id in 1:old_n_elements
            for v in eachvariable(equations)
                old_u[v, .., old_element_id] .= (old_u[v, .., old_element_id] ./
                                                 old_inverse_jacobian[..,
                                                                      old_element_id])
            end
        end

        Trixi.reinitialize_containers!(mesh, equations, dg, cache)

        resize!(u_ode,
                nvariables(equations) * nnodes(dg)^ndims(mesh) * nelements(dg, cache))
        u = Trixi.wrap_array(u_ode, mesh, equations, dg, cache)

        # Loop over all elements in old container and either copy them or coarsen them
        skip = 0
        element_id = 1
        for old_element_id in 1:old_n_elements
            # If skip is non-zero, we just coarsened 2^ndims elements and need to omit the following elements
            if skip > 0
                skip -= 1
                continue
            end

            if to_be_removed[old_element_id]
                # If an element is to be removed, sanity check if the following elements
                # are also marked - otherwise there would be an error in the way the
                # cells/elements are sorted
                @assert all(to_be_removed[old_element_id:(old_element_id + 2^ndims(mesh) - 1)]) "bad cell/element order"

                # Coarsen elements and store solution directly in new data structure
                Trixi.coarsen_elements!(u, element_id, old_u, old_element_id,
                                        adaptor, equations, dg)

                # Before `element_id` is incremented, divide by the new Jacobian and save
                # the result in the parent element
                for v in eachvariable(equations)
                    u[v, .., element_id] .*= (4 .*
                                              cache.elements.inverse_jacobian[..,
                                                                              element_id])
                end

                # Increment `element_id` on the coarsened mesh by one and `skip` = 3 in 2D
                # because 4 children elements become 1 parent element
                element_id += 1
                skip = 2^ndims(mesh) - 1
            else
                # Copy old element data to new element container and remove Jacobian scaling
                for v in eachvariable(equations)
                    u[v, .., element_id] .= (old_u[v, .., old_element_id] .*
                                             old_inverse_jacobian[..,
                                                                  old_element_id])
                end
                # No coarsening occurred, so increment `element_id` on the new mesh by one
                element_id += 1
            end
        end
        # If everything is correct, we should have processed all elements.
        @assert element_id==nelements(dg, cache) + 1 "element_id = $element_id, nelements(dg, cache) = $(nelements(dg, cache))"
    end # GC.@preserve old_u_ode old_inverse_jacobian

    # Apply limiter to ensure admissible solution states
    limiter_shallow_water!(u, equations.threshold_limiter, waterheight,
                           mesh, equations, dg, cache)
    return nothing
end
end # @muladd
