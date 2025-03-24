# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Container data structure (structure-of-arrays style) for DG L2 mortars
# specialized for the shallow water equations. Extra storage is needed to
# store the unprojected parent solution data so that the flux penalty
# can be computed on the mortars directly and then projected back to the parent.
# This ensures that the shallow water solver remains well-balanced on non-conforming meshes.
# For more details on this strategy see Section 3 of the paper
# - Boris Bonev, Jan S. Hesthaven, Francis X. Giraldo, Michal A. Kopera (2018)
#   Discontinuous Galerkin scheme for the spherical shallow water equations
#   with applications to tsunami modeling and prediction
#   [DOI: 10.1016/j.jcp.2018.02.008](https://doi.org/10.1016/j.jcp.2018.02.008)
#
# The positions used in `neighbor_ids` are 1:3 (in 2D) where 1:2 (in 2D)
# are the small elements numbered in z-order and 3 is the large element.
# The solution values on the mortar element are saved in `u`, where `position` is the number
# of the small element that corresponds to the respective part of the mortar element.
# The first dimension `small/large side` and `unprojected` takes 1 for small side,
# 2 for large side, and 3 for the unprojected parent solution.
mutable struct P4estShallowWaterMortarContainer{NDIMS, uEltype <: Real, NDIMSP1, NDIMSP3} <:
               Trixi.AbstractContainer
    u::Array{uEltype, NDIMSP3} # [small/large side, variable, position, i, mortar]
    neighbor_ids::Matrix{Int}             # [position, mortar]
    node_indices::Matrix{NTuple{NDIMS, Symbol}} # [small/large, mortar]

    # internal `resize!`able storage
    _u::Vector{uEltype}
    _neighbor_ids::Vector{Int}
    _node_indices::Vector{NTuple{NDIMS, Symbol}}
end

@inline Trixi.nmortars(mortars::P4estShallowWaterMortarContainer) = size(mortars.neighbor_ids, 2)
@inline Base.ndims(::P4estShallowWaterMortarContainer{NDIMS}) where {NDIMS} = NDIMS

# See explanation of Base.resize! for the element container
function Base.resize!(mortars::P4estShallowWaterMortarContainer, capacity)
    @unpack _u, _neighbor_ids, _node_indices = mortars

    n_dims = ndims(mortars)
    n_nodes = size(mortars.u, 4)
    n_variables = size(mortars.u, 2)

    resize!(_u, 3 * n_variables * 2^(n_dims - 1) * n_nodes^(n_dims - 1) * capacity)
    mortars.u = Trixi.unsafe_wrap(Array, pointer(_u),
                                  (3, n_variables, 2^(n_dims - 1),
                                  ntuple(_ -> n_nodes, n_dims - 1)..., capacity))

    resize!(_neighbor_ids, (2^(n_dims - 1) + 1) * capacity)
    mortars.neighbor_ids = Trixi.unsafe_wrap(Array, pointer(_neighbor_ids),
                                             (2^(n_dims - 1) + 1, capacity))

    resize!(_node_indices, 2 * capacity)
    mortars.node_indices = Trixi.unsafe_wrap(Array, pointer(_node_indices), (2, capacity))

    return nothing
end

# Create mortar container and initialize mortar data.
function Trixi.init_mortars(mesh::Union{P4estMesh, T8codeMesh},
                            equations::ShallowWaterEquationsWetDry2D,
                            basis, elements)
    NDIMS = ndims(elements)
    uEltype = eltype(elements)

    # Initialize container
    n_mortars = Trixi.count_required_surfaces(mesh).mortars

    _u = Vector{uEltype}(undef,
                         3 * nvariables(equations) * 2^(NDIMS - 1) *
                         nnodes(basis)^(NDIMS - 1) * n_mortars)
    u = Trixi.unsafe_wrap(Array, pointer(_u),
                    (3, nvariables(equations), 2^(NDIMS - 1),
                     ntuple(_ -> nnodes(basis), NDIMS - 1)..., n_mortars))

    _neighbor_ids = Vector{Int}(undef, (2^(NDIMS - 1) + 1) * n_mortars)
    neighbor_ids = Trixi.unsafe_wrap(Array, pointer(_neighbor_ids),
                               (2^(NDIMS - 1) + 1, n_mortars))

    _node_indices = Vector{NTuple{NDIMS, Symbol}}(undef, 2 * n_mortars)
    node_indices = Trixi.unsafe_wrap(Array, pointer(_node_indices), (2, n_mortars))

    mortars = P4estShallowWaterMortarContainer{NDIMS, uEltype, NDIMS + 1, NDIMS + 3}(u,
                                                                         neighbor_ids,
                                                                         node_indices,
                                                                         _u,
                                                                         _neighbor_ids,
                                                                         _node_indices)

    if n_mortars > 0
        Trixi.init_mortars!(mortars, mesh)
    end

    return mortars
end
end # @muladd
