# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    VelocityDesingularization()

A stage callback for desingularizing velocities in simulations with small water heights or densities.

This callback is designed to prevent division by very small water heights, which can lead to unphysical
velocities and numerical instabilities. It is applied at each Runge-Kutta stage and modifies the
conservative variables in-place to ensure that the velocity remains well-defined, even in nearly dry
states.

The algorithm works by regularizing the velocity computation and enforcing a minimum water height threshold.
For multi-layer systems, the procedure is applied to each layer individually.
Thresholds to control the desingularization and the minimum water height can be set via 
`threshold_desingularization` and `threshold_limiter` in the equations struct.

Details about the desingularization strategy can be found in Section 2.2 of the paper
- A. Kurganov, G. Petrova (2007)
  A second-order well-balanced positivity preserving central-upwind scheme for the Saint-Venant system
  [doi: 10.4310/CMS.2007.v5.n1.a6](https://dx.doi.org/10.4310/CMS.2007.v5.n1.a6)

The specific implementation is based on Section 3.5 of the paper
- P. Ersing, S. Goldberg, A. Winters (2024)
  Entropy stable hydrostatic reconstruction schemes for shallow water systems
  [doi: 10.1016/j.jcp.2025.113802](https://doi.org/10.1016/j.jcp.2025.113802)
"""
struct VelocityDesingularization end

function (desingularization!::VelocityDesingularization)(u_ode, integrator,
                                                         semi::Trixi.AbstractSemidiscretization,
                                                         t)
    u = Trixi.wrap_array(u_ode, semi)
    Trixi.@trixi_timeit Trixi.timer() "VelocityDesingularization" velocity_desingularization!(u,
                                                                                              Trixi.mesh_equations_solver_cache(semi)...)
end

# Workaround to use this stage callback with time integrators from Trixi.jl
Trixi.init_callback(desingularization!::VelocityDesingularization, semi) = nothing

Trixi.finalize_callback(desingularization!::VelocityDesingularization, semi) = nothing

# Stage callbacks for SimpleSSPRK33 require different arguments (..., stage) instead of (..., semi, t)
function (desingularization!::VelocityDesingularization)(u_ode, integrator, stage)
    semi = integrator.p
    t = integrator.t
    desingularization!(u_ode, integrator, semi, t)
end

@inline function velocity_desingularization!(u, mesh::Trixi.AbstractMesh{1},
                                             equations::ShallowWaterEquations1D,
                                             dg::DGSEM, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)

            h, hv, b = u_node

            # Apply velocity desingularization
            hv = h * (2 * h * hv) /
                 (h^2 + max(h^2, equations.threshold_desingularization))

            if h <= equations.threshold_limiter
                h = equations.threshold_limiter
                hv = zero(eltype(u))
            end

            u_node = SVector(h, hv, b)

            set_node_vars!(u, u_node, equations, dg, i, element)
        end
    end

    return nothing
end

@inline function velocity_desingularization!(u, mesh::Trixi.AbstractMesh{1},
                                             equations::ShallowWaterMultiLayerEquations1D,
                                             dg::DGSEM, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, element)

            h = MVector(waterheight(u_node, equations))
            hv = MVector(momentum(u_node, equations))
            b = u_node[end]

            # Apply velocity desingularization
            hv = h .* (2 * h .* hv) ./
                 (h .^ 2 .+ max.(h .^ 2, equations.threshold_desingularization))

            for i in eachlayer(equations)
                # Ensure positivity and zero velocity at dry states
                if h[i] <= equations.threshold_limiter
                    h[i] = equations.threshold_limiter
                    hv[i] = zero(eltype(u))
                end
            end

            u_node = SVector(h..., hv..., b)

            set_node_vars!(u, u_node, equations, dg, i, element)
        end
    end

    return nothing
end

@inline function velocity_desingularization!(u, mesh::Trixi.AbstractMesh{2},
                                             equations::ShallowWaterEquations2D,
                                             dg::DGSEM, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)

            h, h_v1, h_v2, b = u_node

            # Apply velocity desingularization
            h_v1 = h * (2 * h * h_v1) /
                   (h^2 + max(h^2, equations.threshold_desingularization))
            h_v2 = h * (2 * h * h_v2) /
                   (h^2 + max(h^2, equations.threshold_desingularization))

            if h <= equations.threshold_limiter
                h = equations.threshold_limiter
                h_v1 = zero(eltype(u))
                h_v2 = zero(eltype(u))
            end

            u_node = SVector(h, h_v1, h_v2, b)

            set_node_vars!(u, u_node, equations, dg, i, j, element)
        end
    end

    return nothing
end

@inline function velocity_desingularization!(u, mesh::Trixi.AbstractMesh{2},
                                             equations::ShallowWaterMultiLayerEquations2D,
                                             dg::DGSEM, cache)
    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)

            h = MVector(waterheight(u_node, equations))
            h_v1 = MVector(momentum(u_node, equations)[1])
            h_v2 = MVector(momentum(u_node, equations)[2])
            b = u_node[end]

            # Apply velocity desingularization
            h_v1 = h .* (2 * h .* h_v1) ./
                   (h .^ 2 .+ max.(h .^ 2, equations.threshold_desingularization))
            h_v2 = h .* (2 * h .* h_v2) ./
                   (h .^ 2 .+ max.(h .^ 2, equations.threshold_desingularization))

            for i in eachlayer(equations)
                # Ensure positivity and zero velocity at dry states
                if h[i] <= equations.threshold_limiter
                    h[i] = equations.threshold_limiter
                    h_v1[i] = zero(eltype(u))
                    h_v2[i] = zero(eltype(u))
                end
            end

            u_node = SVector(h..., h_v1..., h_v2..., b)

            set_node_vars!(u, u_node, equations, dg, i, j, element)
        end
    end

    return nothing
end
end # @muladd
