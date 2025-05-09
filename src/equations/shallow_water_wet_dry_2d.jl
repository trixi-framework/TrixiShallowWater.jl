# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent
@doc raw"""
    ShallowWaterEquations2D(; gravity, H0 = 0, threshold_limiter = nothing, threshold_wet = nothing, 
                              threshold_partially_wet = nothing, threshold_desingularization = nothing)

Shallow water equations (SWE) in two space dimensions. The equations are given by
```math
\begin{aligned}
  \frac{\partial h}{\partial t} + \frac{\partial}{\partial x}(h v_1)
    + \frac{\partial}{\partial y}(h v_2) &= 0 \\
    \frac{\partial}{\partial t}(h v_1) + \frac{\partial}{\partial x}\left(h v_1^2 + \frac{g}{2}h^2\right)
    + \frac{\partial}{\partial y}(h v_1 v_2) + g h \frac{\partial b}{\partial x} &= 0 \\
    \frac{\partial}{\partial t}(h v_2) + \frac{\partial}{\partial x}(h v_1 v_2)
    + \frac{\partial}{\partial y}\left(h v_2^2 + \frac{g}{2}h^2\right) + g h \frac{\partial b}{\partial y} &= 0.
\end{aligned}
```
The unknown quantities of the SWE are the water height ``h`` and the velocities ``\mathbf{v} = (v_1, v_2)^T``.
The gravitational acceleration is denoted by `g` and the (possibly) variable bottom topography function ``b(x,y)``.
Conservative variable water height ``h`` is measured from the bottom topography ``b``, therefore one
also defines the total water height as ``H = h + b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water height that
is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

Also, there are four thresholds which prevent numerical problems as well as instabilities. None of them
have to be passed, as default values are defined within the struct. The first one, `threshold_limiter`, is
used in [`PositivityPreservingLimiterShallowWater`](@ref) on the water height, as a (small) shift on the initial
condition and cutoff before the next time step. The second one, `threshold_wet`, is applied on the water height to
define when the flow is "wet" before calculating the numerical flux. A third
`threshold_partially_wet` is applied on the water height to define "partially wet" elements in
[`IndicatorHennemannGassnerShallowWater`](@ref), that are then calculated with a pure FV method to
ensure well-balancedness. Lastly, `threshold_desingularization` is used in [`PositivityPreservingLimiterShallowWater`](@ref)
for the velocity desingularization procedure. For `Float64` no threshold needs to be passed, as default values are
defined within the struct. For other number formats  `threshold_partially_wet`
and `threshold_desingularization` must be provided.

The bottom topography function ``b(x,y)`` is set inside the initial condition routine
for a particular problem setup. To test the conservative form of the SWE one can set the bottom topography
variable `b` to zero.

In addition to the unknowns, TrixiShallowWater.jl currently stores the bottom topography values at the approximation points
despite being fixed in time. This is done for convenience of computing the bottom topography gradients
on the fly during the approximation as well as computing auxiliary quantities like the total water height ``H``
or the entropy variables.
This affects the implementation and use of these equations in various ways:
* The flux values corresponding to the bottom topography must be zero.
* The bottom topography values must be included when defining initial conditions, boundary conditions or
  source terms.
* [`Trixi.AnalysisCallback`](@extref) analyzes this variable.
* Trixi.jl's visualization tools will visualize the bottom topography by default.

References for the SWE are many but a good introduction is available in Chapter 13 of the book:
- Randall J. LeVeque (2002)
  Finite Volume Methods for Hyperbolic Problems
  [DOI: 10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253)
"""
struct ShallowWaterEquations2D{RealT <: Real} <:
       Trixi.AbstractShallowWaterEquations{2, 4}
    gravity::RealT # gravitational acceleration
    H0::RealT      # constant "lake-at-rest" total water height
    # `threshold_limiter` used in `PositivityPreservingLimiterShallowWater` on water height,
    # as a (small) shift on the initial condition and cutoff before the next time step.
    # Default is 500*eps() which in double precision is ≈1e-13.
    threshold_limiter::RealT
    # `threshold_wet` applied on water height to define when the flow is "wet"
    # before calculating the numerical flux.
    # Default is 5*eps() which in double precision is ≈1e-15.
    threshold_wet::RealT
    # `threshold_partially_wet` used in `IndicatorHennemannGassnerShallowWater` on the water height
    # to define "partially wet" elements. Those elements are calculated with a pure FV method to
    # ensure well-balancedness. Default in double precision is 1e-4.
    threshold_partially_wet::RealT
    # `threshold_desingularization` used in the velocity desingularization procedure, to avoid
    # division by small numbers. Default in double precision is 1e-10.
    threshold_desingularization::RealT
end

# Allow for flexibility to set the gravitational acceleration within an elixir depending on the
# application where `gravity=1.0` or `gravity=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases.
# Strict default values for thresholds that performed well in many numerical experiments
function ShallowWaterEquations2D(; gravity, H0 = zero(gravity),
                                 threshold_limiter = nothing,
                                 threshold_wet = nothing,
                                 threshold_partially_wet = nothing,
                                 threshold_desingularization = nothing)
    T = promote_type(typeof(gravity), typeof(H0))
    if threshold_limiter === nothing
        threshold_limiter = 500 * eps(T)
    end
    if threshold_wet === nothing
        threshold_wet = 5 * eps(T)
    end
    if threshold_partially_wet === nothing
        threshold_partially_wet = default_threshold_partially_wet(T)
    end
    if threshold_desingularization === nothing
        threshold_desingularization = default_threshold_desingularization(T)
    end

    ShallowWaterEquations2D(gravity, H0, threshold_limiter,
                            threshold_wet, threshold_partially_wet,
                            threshold_desingularization)
end

Trixi.have_nonconservative_terms(::ShallowWaterEquations2D) = True()

Trixi.varnames(::typeof(cons2cons), ::ShallowWaterEquations2D) = ("h", "h_v1",
                                                                  "h_v2", "b")
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
Trixi.varnames(::typeof(cons2prim), ::ShallowWaterEquations2D) = ("H", "v1", "v2",
                                                                  "b")

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterEquations2D)

A smooth initial condition used for convergence tests in combination with
[`source_terms_convergence_test`](@ref)
(and [`Trixi.BoundaryConditionDirichlet`](@extref) in non-periodic domains).
"""
function Trixi.initial_condition_convergence_test(x, t,
                                                  equations::ShallowWaterEquations2D)
    # some constants are chosen such that the function is periodic on the domain [0,sqrt(2)]^2
    RealT = eltype(x)
    c = 7
    omega_x = 2 * convert(RealT, pi) * sqrt(convert(RealT, 2))
    omega_t = 2 * convert(RealT, pi)

    x1, x2 = x

    H = c + cos(omega_x * x1) * sin(omega_x * x2) * cos(omega_t * t)
    v1 = 0.5f0
    v2 = 1.5f0
    b = 2 + 0.5f0 * sinpi(sqrt(convert(RealT, 2)) * x1) +
        0.5f0 * sinpi(sqrt(convert(RealT, 2)) * x2)
    return prim2cons(SVector(H, v1, v2, b), equations)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterEquations2D)

Source terms used for convergence tests in combination with
[`initial_condition_convergence_test`](@ref)
(and [`Trixi.BoundaryConditionDirichlet`](@extref) in non-periodic domains).

This manufactured solution source term is specifically designed for the bottom topography function
`b(x,y) = 2 + 0.5 * sin(sqrt(2)*pi*x) + 0.5 * sin(sqrt(2)*pi*y)`
as defined in [`initial_condition_convergence_test`](@ref).
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterEquations2D)
    # Same settings as in `initial_condition_convergence_test`. Some derivatives simplify because
    # for this manufactured solution velocities are taken to be constants
    RealT = eltype(u)
    c = 7
    omega_x = 2 * convert(RealT, pi) * sqrt(convert(RealT, 2))
    omega_t = 2 * convert(RealT, pi)
    omega_b = sqrt(convert(RealT, 2)) * convert(RealT, pi)
    v1 = 0.5f0
    v2 = 1.5f0

    x1, x2 = x

    sinX, cosX = sincos(omega_x * x1)
    sinY, cosY = sincos(omega_x * x2)
    sinT, cosT = sincos(omega_t * t)

    H = c + cosX * sinY * cosT
    H_x = -omega_x * sinX * sinY * cosT
    H_y = omega_x * cosX * cosY * cosT
    # this time derivative for the water height exploits that the bottom topography is
    # fixed in time such that H_t = (h+b)_t = h_t + 0
    H_t = -omega_t * cosX * sinY * sinT

    # bottom topography and its gradient
    b = 2 + 0.5f0 * sinpi(sqrt(convert(RealT, 2)) * x1) +
        0.5f0 * sinpi(sqrt(convert(RealT, 2)) * x2)
    tmp1 = 0.5f0 * omega_b
    b_x = tmp1 * cos(omega_b * x1)
    b_y = tmp1 * cos(omega_b * x2)

    du1 = H_t + v1 * (H_x - b_x) + v2 * (H_y - b_y)
    du2 = v1 * du1 + equations.gravity * (H - b) * H_x
    du3 = v2 * du1 + equations.gravity * (H - b) * H_y
    return SVector(du1, du2, du3, 0)
end

"""
    initial_condition_weak_blast_wave(x, t, equations::ShallowWaterEquations2D)

A weak blast wave discontinuity useful for testing, e.g., total energy conservation.
Note for the shallow water equations to the total energy acts as a mathematical entropy function.
"""
function Trixi.initial_condition_weak_blast_wave(x, t,
                                                 equations::ShallowWaterEquations2D)
    # Set up polar coordinates
    RealT = eltype(x)
    inicenter = SVector(convert(RealT, 0.7), convert(RealT, 0.7))
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    phi = atan(y_norm, x_norm)
    sin_phi, cos_phi = sincos(phi)

    # Calculate primitive variables
    H = r > 0.5f0 ? 3.25f0 : 4.0f0
    v1 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * cos_phi
    v2 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * sin_phi
    b = 0 # by default assume there is no bottom topography

    return prim2cons(SVector(H, v1, v2, b), equations)
end

"""
    boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,
                                 equations::ShallowWaterEquations2D)
Create a boundary state by reflecting the normal velocity component and keep
the tangential velocity component unchanged. The boundary water height is taken from
the internal value.
For details see Section 9.2.5 of the book:
- Eleuterio F. Toro (2001)
  Shock-Capturing Methods for Free-Surface Shallow Flows
  1st edition
  ISBN 0471987662
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner,
                                                    normal_direction::AbstractVector,
                                                    x, t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

"""
    boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                 surface_flux_function, equations::ShallowWaterEquations2D)

Should be used together with [`Trixi.TreeMesh`](@extref).
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation,
                                                    direction, x, t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])
    else # orientation == 2
        u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4])
    end

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
                                                     equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
        noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
                                                     equations)
    end

    return flux, noncons_flux
end

"""
    BoundaryConditionWaterHeight(h_boundary, equations::ShallowWaterEquations2D)

Create a boundary condition that uses `h_boundary` to specify a fixed water height at the
boundary and extrapolates the velocity from the incoming Riemann invariant.

The external water height `h_boundary` can be specified as a constant value or as a function of time, e.g.
```julia
   BoundaryConditionWaterHeight(h_boundary, equations))
   BoundaryConditionWaterHeight(t -> h_boundary(t), equations))
```

More details can be found in the paper:
- Lixiang Song, Jianzhong Zhou, Jun Guo, Qiang Zou, Yi Liu (2011)
  A robust well-balanced finite volume model for shallow water flows
  with wetting and drying over irregular terrain
  [doi: 10.1016/j.advwatres.2011.04.017](https://doi.org/10.1016/j.advwatres.2011.04.017)

!!! warning "Experimental code"
    This is an experimental feature and can change any time.
"""
function BoundaryConditionWaterHeight(h_boundary::Real,
                                      equations::ShallowWaterEquations2D{RealT}) where {RealT}
    # Convert function output to the correct type
    h_boundary = convert(RealT, h_boundary)
    return BoundaryConditionWaterHeight(t -> h_boundary)
end

function BoundaryConditionWaterHeight(h_boundary::Function,
                                      equations::ShallowWaterEquations2D{RealT}) where {RealT}
    # Check if the function output is of the correct type
    if !(typeof(h_boundary(one(RealT))) == RealT)
        throw(ArgumentError("Boundary value functions must return a value of type $(RealT)"))
    end
    return BoundaryConditionWaterHeight(t -> h_boundary(t))
end

# Version for `TreeMesh`
function (boundary_condition::BoundaryConditionWaterHeight)(u_inner,
                                                            orientation,
                                                            direction, x, t,
                                                            surface_flux_functions,
                                                            equations::ShallowWaterEquations2D)
    # Extract the gravitational acceleration
    g = equations.gravity

    # Get the water height and velocity from the inner state
    h_inner = waterheight(u_inner, equations)
    v1_inner, v2_inner = velocity(u_inner, equations)

    # Extract the external water height from the boundary condition
    h_boundary = boundary_condition.h_boundary(t)

    # Calculate the boundary state based on the direction.
    # To extrapolate the external velocity assume that the Riemann invariant remains constant across
    # the incoming characteristic. In the case of inflow we assume that the tangential velocity at
    # the boundary is zero.
    if direction == 1 # x-
        v1_boundary = v1_inner +
                      2 * (sqrt(g * h_boundary) - sqrt(g * h_inner))
        v1_boundary > 0 ? hv2_boundary = zero(u_inner[3]) : hv2_boundary = u_inner[3]
        u_boundary = SVector(h_boundary, h_boundary * v1_boundary, hv2_boundary,
                             u_inner[4])
    elseif direction == 2 # x+
        v1_boundary = v1_inner -
                      2 * (sqrt(g * h_boundary) - sqrt(g * h_inner))
        v1_boundary < 0 ? hv2_boundary = zero(u_inner[3]) : hv2_boundary = u_inner[3]
        u_boundary = SVector(h_boundary, h_boundary * v1_boundary, hv2_boundary,
                             u_inner[4])
    elseif direction == 3 # y-
        v2_boundary = v2_inner +
                      2 * (sqrt(g * h_boundary) - sqrt(g * h_inner))
        v2_boundary > 0 ? hv1_boundary = zero(u_inner[2]) : hv1_boundary = u_inner[2]
        u_boundary = SVector(h_boundary, hv1_boundary, h_boundary * v2_boundary,
                             u_inner[4])
    elseif direction == 4 # y+
        v2_boundary = v2_inner -
                      2 * (sqrt(g * h_boundary) - sqrt(g * h_inner))
        v2_boundary < 0 ? hv1_boundary = zero(u_inner[2]) : hv1_boundary = u_inner[2]
        u_boundary = SVector(h_boundary, hv1_boundary, h_boundary * v2_boundary,
                             u_inner[4])
    end

    # Evaluate the conservative flux at the boundary
    flux = Trixi.flux(u_boundary, orientation, equations)

    # Return the conservative and nonconservative fluxes.
    # The nonconservative part is zero as we assume a constant bottom topography at the boundary.
    return (flux, zero(u_inner))
end

# Version for `UnstructuredMesh2D` and `P4estMesh`
function (boundary_condition::BoundaryConditionWaterHeight)(u_inner,
                                                            normal_direction,
                                                            x, t,
                                                            surface_flux_functions,
                                                            equations::ShallowWaterEquations2D)
    # Extract the gravitational acceleration
    g = equations.gravity

    # Normalized normal vector
    normal = normal_direction / Trixi.norm(normal_direction)

    # Apply the rotation that maps `normal` onto the x-axis to `u_inner`.
    u_rotated = Trixi.rotate_to_x(u_inner, normal, equations)

    # Get the water height and velocity from the inner state
    h_inner = waterheight(u_rotated, equations)
    v_inner_normal, _ = velocity(u_rotated, equations)

    # Extract the external water height from the boundary condition
    h_boundary = boundary_condition.h_boundary(t)

    # Calculate the boundary state in the rotated coordinate system.
    # To extrapolate the external velocity assume that the Riemann invariant remains constant across
    # the incoming characteristic. In the case of inflow we assume that the tangential velocity at
    # the boundary is zero.
    v_boundary_normal = v_inner_normal - 2 * (sqrt(g * h_boundary) - sqrt(g * h_inner))
    v_boundary_normal < 0 ? hv_boundary_tangential = zero(u_rotated[3]) :
    hv_boundary_tangential = u_rotated[3]

    u_boundary = SVector(h_boundary, h_boundary * v_boundary_normal,
                         hv_boundary_tangential, u_rotated[4])

    # Compute the boundary flux in the rotated coordinate system.
    flux = Trixi.flux(u_boundary, 1, equations)

    # Apply the back-rotation that maps the x-axis onto `normal` to the boundary flux.
    flux = Trixi.rotate_from_x(flux, normal, equations) * Trixi.norm(normal_direction)

    # Return the conservative and nonconservative fluxes. The nonconservative part is zero as we assume
    # a constant bottom topography at the boundary.
    return (flux, zero(u_inner))
end

"""
    BoundaryConditionMomentum(hv1_boundary, hv2_boundary, equations::ShallowWaterEquations2D)

Create a boundary condition that sets a fixed momentum in x- and y- directions, `hv1_boundary`
and `hv2_boundary`, at the boundary and extrapolates the water height `h_boundary` from the incoming
Riemann invariant.

The external momentum can be specified as a constant value or as a function of time, e.g.
```julia
   BoundaryConditionMomentum(hv1_boundary, hv2_boundary, equations)
   BoundaryConditionMomentum(t -> hv1_boundary(t), t -> hv2_boundary(t), equations)
```

More details can be found in the paper:
- Lixiang Song, Jianzhong Zhou, Jun Guo, Qiang Zou, Yi Liu (2011)
  A robust well-balanced finite volume model for shallow water flows
  with wetting and drying over irregular terrain
  [doi: 10.1016/j.advwatres.2011.04.017](https://doi.org/10.1016/j.advwatres.2011.04.017)

!!! warning "Experimental code"
    This is an experimental feature and can change any time.
"""
function BoundaryConditionMomentum(hv1_boundary::Real, hv2_boundary::Real,
                                   equations::ShallowWaterEquations2D{RealT}) where {RealT}
    # Convert function output to the correct type
    hv1_boundary = convert(RealT, hv1_boundary)
    hv2_boundary = convert(RealT, hv2_boundary)
    return BoundaryConditionMomentum((t -> (hv1_boundary, hv2_boundary)))
end

function BoundaryConditionMomentum(hv1_boundary::Function, hv2_boundary::Function,
                                   equations::ShallowWaterEquations2D{RealT}) where {RealT}
    # Check if the function output is of the correct type
    if !(typeof(hv1_boundary(one(RealT))) == RealT &&
         typeof(hv2_boundary(one(RealT))) == RealT)
        throw(ArgumentError("Boundary value functions must return a value of type $(RealT)"))
    end
    return BoundaryConditionMomentum(t -> (hv1_boundary(t), hv2_boundary(t)))
end

# Version for `TreeMesh`
function (boundary_condition::BoundaryConditionMomentum)(u_inner,
                                                         orientation,
                                                         direction, x, t,
                                                         surface_flux_functions,
                                                         equations::ShallowWaterEquations2D)
    # Extract the gravitational acceleration
    g = equations.gravity

    # Get the water height and velocity from the inner state
    h_inner = waterheight(u_inner, equations)
    v1_inner, v2_inner = velocity(u_inner, equations)

    # Extract the external momentum from the boundary condition
    hv1_boundary, hv2_boundary = boundary_condition.hv_boundary(t)

    # Calculate the boundary state based on the direction.
    # To extrapolate the external water height `h_boundary` assume that the Riemann invariant remains
    # constant across the incoming characteristic.
    # Requires one to solve for the roots of a nonlinear function, see Eq. (52) in the reference above.
    # For convenience we substitute x = h_boundary and solve for x, using the Steffensen method.
    if direction == 1 # x-
        fx = ZeroProblem(x -> 2 * sqrt(g) * x^(3 / 2) +
                              (v1_inner - 2 * sqrt(g * h_inner)) * x - hv1_boundary,
                         h_inner)
        h_boundary = solve(fx, Order2())
        if hv1_boundary > 0
            u_boundary = SVector(h_boundary, hv1_boundary, hv2_boundary, u_inner[4])
        else
            u_boundary = SVector(h_boundary, hv1_boundary, u_inner[3], u_inner[4])
        end
    elseif direction == 2 # x+
        fx = ZeroProblem(x -> 2 * sqrt(g) * x^(3 / 2) -
                              (v1_inner + 2 * sqrt(g * h_inner)) * x + hv1_boundary,
                         h_inner)
        h_boundary = solve(fx, Order2())
        if hv1_boundary < 0
            u_boundary = SVector(h_boundary, hv1_boundary, hv2_boundary, u_inner[4])
        else
            u_boundary = SVector(h_boundary, hv1_boundary, u_inner[3], u_inner[4])
        end
    elseif direction == 3 # y-
        fx = ZeroProblem(x -> 2 * sqrt(g) * x^(3 / 2) +
                              (v2_inner - 2 * sqrt(g * h_inner)) * x - hv2_boundary,
                         h_inner)
        h_boundary = solve(fx, Order2())
        if hv2_boundary > 0
            u_boundary = SVector(h_boundary, hv1_boundary, hv2_boundary, u_inner[4])
        else
            u_boundary = SVector(h_boundary, u_inner[2], hv2_boundary, u_inner[4])
        end
    elseif direction == 4 # y+
        fx = ZeroProblem(x -> 2 * sqrt(g) * x^(3 / 2) -
                              (v2_inner + 2 * sqrt(g * h_inner)) * x + hv2_boundary,
                         h_inner)
        h_boundary = solve(fx, Order2())
        if hv2_boundary < 0
            u_boundary = SVector(h_boundary, hv1_boundary, hv2_boundary, u_inner[4])
        else
            u_boundary = SVector(h_boundary, u_inner[2], hv2_boundary, u_inner[4])
        end
    end

    # Evaluate the conservative flux at the boundary
    flux = Trixi.flux(u_boundary, orientation, equations)

    # Return the conservative and nonconservative fluxes.
    # The nonconservative part is zero as we assume a constant bottom topography at the boundary.
    return (flux, zero(u_inner))
end

# Version for `UnstructuredMesh2D` and `P4estMesh`
function (boundary_condition::BoundaryConditionMomentum)(u_inner,
                                                         normal_direction,
                                                         x, t,
                                                         surface_flux_functions,
                                                         equations::ShallowWaterEquations2D)

    # Extract the gravitational acceleration
    g = equations.gravity

    # Normalized normal vector
    normal = normal_direction / Trixi.norm(normal_direction)

    # Apply the rotation that maps `normal` onto the x-axis to `u_inner` and `hv_boundary``.
    u_rotated = Trixi.rotate_to_x(u_inner, normal, equations)

    # Get the water height and velocity from the inner state
    h_inner = waterheight(u_rotated, equations)
    v_inner_normal, _ = velocity(u_rotated, equations)

    # Extract the external momentum from the boundary condition
    hv1_boundary, hv2_boundary = boundary_condition.hv_boundary(t)

    hv_boundary_normal = hv1_boundary * normal[1] + hv2_boundary * normal[2]
    hv_boundary_tangential = -hv1_boundary * normal[2] + hv2_boundary * normal[1]

    # Calculate the boundary state in the rotated coordinate system.
    # To extrapolate the external water height `h_boundary` assume that the Riemann invariant remains
    # constant across the incoming characteristic.
    # Requires one to solve for the roots of a nonlinear function, see Eq. (52) in the reference above.
    # For convenience we substitute x = h_boundary and solve for x.
    fx = ZeroProblem(x -> 2 * sqrt(g) * x^(3 / 2) -
                          (v_inner_normal + 2 * sqrt(g * h_inner)) * x +
                          hv_boundary_normal, h_inner)
    h_boundary = solve(fx, Order2())

    hv_boundary_normal < 0 ? nothing : hv_boundary_tangential = u_rotated[3]

    u_boundary = SVector(h_boundary, hv_boundary_normal, hv_boundary_tangential,
                         u_inner[4])

    # Compute the boundary flux in the rotated coordinate system.
    flux = Trixi.flux(u_boundary, 1, equations)

    # Apply the back-rotation that maps the x-axis onto `normal` to the boundary flux.
    flux = Trixi.rotate_from_x(flux, normal, equations) * Trixi.norm(normal_direction)

    # Return the conservative and nonconservative fluxes. The nonconservative part is zero as we assume
    # a constant bottom topography at the boundary.
    return (flux, zero(u_inner))
end

# Calculate 1D flux for a single point
# Note, the bottom topography has no flux
@inline function Trixi.flux(u, orientation::Integer,
                            equations::ShallowWaterEquations2D)
    h, h_v1, h_v2, _ = u
    v1, v2 = velocity(u, equations)

    p = 0.5f0 * equations.gravity * h^2
    if orientation == 1
        f1 = h_v1
        f2 = h_v1 * v1 + p
        f3 = h_v1 * v2
    else
        f1 = h_v2
        f2 = h_v2 * v1
        f3 = h_v2 * v2 + p
    end
    return SVector(f1, f2, f3, 0)
end

# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized and the bottom topography has no flux
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::ShallowWaterEquations2D)
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    h_v_normal = h * v_normal
    p = 0.5f0 * equations.gravity * h^2

    f1 = h_v_normal
    f2 = h_v_normal * v1 + p * normal_direction[1]
    f3 = h_v_normal * v2 + p * normal_direction[2]
    return SVector(f1, f2, f3, 0)
end

"""
    flux_nonconservative_wintermeyer_etal(u_ll, u_rr, orientation::Integer,
                                          equations::ShallowWaterEquations2D)
    flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                          normal_direction::AbstractVector,
                                          equations::ShallowWaterEquations2D)

Non-symmetric two-point volume flux discretizing the nonconservative (source) term
that contains the gradient of the bottom topography [`ShallowWaterEquations2D`](@ref).

For the `surface_flux` either [`flux_wintermeyer_etal`](@ref) or [`flux_fjordholm_etal`](@ref) can
be used to ensure well-balancedness and entropy conservation.

Further details are available in the papers:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on
  curvilinear meshes
  [DOI: 10.48550/arXiv.2306.12699](https://doi.org/10.48550/arXiv.2306.12699)
"""
@inline function Trixi.flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                             orientation::Integer,
                                                             equations::ShallowWaterEquations2D)
    # Pull the necessary left and right state information
    h_ll = waterheight(u_ll, equations)
    b_jump = u_rr[4] - u_ll[4]

    # Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
    if orientation == 1
        f = SVector(0, equations.gravity * h_ll * b_jump, 0, 0)
    else # orientation == 2
        f = SVector(0, 0, equations.gravity * h_ll * b_jump, 0)
    end
    return f
end

@inline function Trixi.flux_nonconservative_wintermeyer_etal(u_ll, u_rr,
                                                             normal_direction::AbstractVector,
                                                             equations::ShallowWaterEquations2D)
    # Pull the necessary left and right state information
    h_ll = waterheight(u_ll, equations)
    b_jump = u_rr[4] - u_ll[4]

    # Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
    return SVector(0,
                   normal_direction[1] * equations.gravity * h_ll * b_jump,
                   normal_direction[2] * equations.gravity * h_ll * b_jump,
                   0)
end

"""
    flux_nonconservative_fjordholm_etal(u_ll, u_rr, orientation::Integer,
                                        equations::ShallowWaterEquations2D)
    flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                        normal_direction::AbstractVector,
                                        equations::ShallowWaterEquations2D)

Non-symmetric two-point surface flux discretizing the nonconservative (source) term of
that contains the gradient of the bottom topography [`ShallowWaterEquations2D`](@ref).

This flux can be used together with [`flux_fjordholm_etal`](@ref) at interfaces to ensure entropy
conservation and well-balancedness.

Further details for the original finite volume formulation are available in
- Ulrik S. Fjordholm, Siddhartha Mishra and Eitan Tadmor (2011)
  Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography
  [DOI: 10.1016/j.jcp.2011.03.042](https://doi.org/10.1016/j.jcp.2011.03.042)
and for curvilinear 2D case in the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
"""
@inline function Trixi.flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                           orientation::Integer,
                                                           equations::ShallowWaterEquations2D)
    # Pull the necessary left and right state information
    h_ll, _, _, b_ll = u_ll
    h_rr, _, _, b_rr = u_rr

    h_average = 0.5f0 * (h_ll + h_rr)
    b_jump = b_rr - b_ll

    # Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
    if orientation == 1
        f = SVector(0,
                    equations.gravity * h_average * b_jump,
                    0, 0)
    else # orientation == 2
        f = SVector(0, 0,
                    equations.gravity * h_average * b_jump,
                    0)
    end

    return f
end

@inline function Trixi.flux_nonconservative_fjordholm_etal(u_ll, u_rr,
                                                           normal_direction::AbstractVector,
                                                           equations::ShallowWaterEquations2D)
    # Pull the necessary left and right state information
    h_ll, _, _, b_ll = u_ll
    h_rr, _, _, b_rr = u_rr

    h_average = 0.5f0 * (h_ll + h_rr)
    b_jump = b_rr - b_ll

    # Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
    f2 = normal_direction[1] * equations.gravity * h_average * b_jump
    f3 = normal_direction[2] * equations.gravity * h_average * b_jump

    # First and last equations do not have a nonconservative flux
    f1 = f4 = 0

    return SVector(f1, f2, f3, f4)
end

"""
    hydrostatic_reconstruction_audusse_etal(u_ll, u_rr, orientation_or_normal_direction,
                                            equations::ShallowWaterEquations2D)

A particular type of hydrostatic reconstruction on the water height to guarantee well-balancedness
for a general bottom topography [`ShallowWaterEquations2D`](@ref). The reconstructed solution states
`u_ll_star` and `u_rr_star` variables are used to evaluate the surface numerical flux at the interface.
Use in combination with the generic numerical flux routine [`Trixi.FluxHydrostaticReconstruction`](@extref).

Further details for the hydrostatic reconstruction and its motivation can be found in
- Emmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004)
  A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
  [DOI: 10.1137/S1064827503431090](https://doi.org/10.1137/S1064827503431090)
"""
@inline function hydrostatic_reconstruction_audusse_etal(u_ll, u_rr,
                                                         equations::ShallowWaterEquations2D)
    # Unpack left and right water heights and bottom topographies
    h_ll, _, _, b_ll = u_ll
    h_rr, _, _, b_rr = u_rr

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Compute the reconstructed water heights
    h_ll_star = max(0, h_ll + b_ll - max(b_ll, b_rr))
    h_rr_star = max(0, h_rr + b_rr - max(b_ll, b_rr))

    # Create the conservative variables using the reconstruted water heights
    u_ll_star = SVector(h_ll_star, h_ll_star * v1_ll, h_ll_star * v2_ll, b_ll)
    u_rr_star = SVector(h_rr_star, h_rr_star * v1_rr, h_rr_star * v2_rr, b_rr)

    return u_ll_star, u_rr_star
end
# TODO: Adapt this to wetting and drying

"""
    hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquations2D)

A particular type of hydrostatic reconstruction of the water height to guarantee well-balancedness
for a general bottom topography of the [`ShallowWaterEquations2D`](@ref). The reconstructed solution states
`u_ll_star` and `u_rr_star` variables are then used to evaluate the surface numerical flux at the interface.
The key idea is a linear reconstruction of the bottom and water height at the interfaces using subcells.
Use in combination with the generic numerical flux routine [`Trixi.FluxHydrostaticReconstruction`](@extref).

Further details on this hydrostatic reconstruction and its motivation can be found in
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI:10.1137/15M1053074](https://dx.doi.org/10.1137/15M1053074)
"""
@inline function hydrostatic_reconstruction_chen_noelle(u_ll, u_rr,
                                                        equations::ShallowWaterEquations2D)
    # Unpack left and right water heights and bottom topographies
    h_ll, _, _, b_ll = u_ll
    h_rr, _, _, b_rr = u_rr

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    H_ll = b_ll + h_ll
    H_rr = b_rr + h_rr

    b_star = min(max(b_ll, b_rr), min(H_ll, H_rr))

    # Compute the reconstructed water heights
    h_ll_star = min(H_ll - b_star, h_ll)
    h_rr_star = min(H_rr - b_star, h_rr)

    # Set the water height to be at least the value stored in the variable threshold after
    # the hydrostatic reconstruction is applied and before the numerical flux is calculated
    # to avoid numerical problem with arbitrary small values. Interfaces with a water height
    # lower or equal to the threshold can be declared as dry.
    # The default value for `threshold_wet` is ≈5*eps(), or 1e-15 in double precision, is set
    # in the `ShallowWaterEquations2D` struct. This threshold value can be changed in the constructor
    # call of this equation struct in an elixir.
    threshold = equations.threshold_wet

    if (h_ll_star <= threshold)
        h_ll_star = threshold
        v1_ll = zero(v1_ll)
        v2_ll = zero(v2_ll)
    end

    if (h_rr_star <= threshold)
        h_rr_star = threshold
        v1_rr = zero(v1_rr)
        v2_rr = zero(v2_rr)
    end

    # Create the conservative variables using the reconstruted water heights
    u_ll_star = SVector(h_ll_star, h_ll_star * v1_ll, h_ll_star * v2_ll, b_ll)
    u_rr_star = SVector(h_rr_star, h_rr_star * v1_rr, h_rr_star * v2_rr, b_rr)

    return u_ll_star, u_rr_star
end

"""
    flux_nonconservative_audusse_etal(u_ll, u_rr, orientation::Integer,
                                      equations::ShallowWaterEquations2D)
    flux_nonconservative_audusse_etal(u_ll, u_rr,
                                      normal_direction::AbstractVector,
                                      equations::ShallowWaterEquations2D)

Non-symmetric two-point surface flux that discretizes the nonconservative (source) term.
The discretization uses the [`hydrostatic_reconstruction_audusse_etal`](@ref) on the conservative
variables.

This hydrostatic reconstruction ensures that the finite volume numerical fluxes remain
well-balanced for discontinuous bottom topographies [`ShallowWaterEquations2D`](@ref).
Should be used together with [`Trixi.FluxHydrostaticReconstruction`](@extref) and
[`hydrostatic_reconstruction_audusse_etal`](@ref) in the surface flux to ensure consistency.

Further details for the hydrostatic reconstruction and its motivation can be found in
- Emmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004)
  A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
  [DOI: 10.1137/S1064827503431090](https://doi.org/10.1137/S1064827503431090)
"""
@inline function flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                   orientation::Integer,
                                                   equations::ShallowWaterEquations2D)
    # Pull the water height and bottom topography on the left
    h_ll, _, _, b_ll = u_ll

    # Create the hydrostatic reconstruction for the left solution state
    u_ll_star, _ = hydrostatic_reconstruction_audusse_etal(u_ll, u_rr, equations)

    # Copy the reconstructed water height for easier to read code
    h_ll_star = u_ll_star[1]

    if orientation == 1
        f = SVector(0,
                    equations.gravity * (h_ll^2 - h_ll_star^2),
                    0, 0)
    else # orientation == 2
        f = SVector(0, 0,
                    equations.gravity * (h_ll^2 - h_ll_star^2),
                    0)
    end

    return f
end

@inline function flux_nonconservative_audusse_etal(u_ll, u_rr,
                                                   normal_direction::AbstractVector,
                                                   equations::ShallowWaterEquations2D)
    # Pull the water height and bottom topography on the left
    h_ll, _, _, b_ll = u_ll

    # Create the hydrostatic reconstruction for the left solution state
    u_ll_star, _ = hydrostatic_reconstruction_audusse_etal(u_ll, u_rr, equations)

    # Copy the reconstructed water height for easier to read code
    h_ll_star = u_ll_star[1]

    f2 = normal_direction[1] * equations.gravity * (h_ll^2 - h_ll_star^2)
    f3 = normal_direction[2] * equations.gravity * (h_ll^2 - h_ll_star^2)

    # First and last equations do not have a nonconservative flux
    f1 = f4 = 0

    return SVector(f1, f2, f3, f4)
end

"""
    flux_nonconservative_chen_noelle(u_ll, u_rr,
                                     orientation::Integer,
                                     equations::ShallowWaterEquations2D)
    flux_nonconservative_chen_noelle(u_ll, u_rr,
                                     normal_direction::AbstractVector,
                                     equations::ShallowWaterEquations2D)

Non-symmetric two-point surface flux that discretizes the nonconservative (source) term.
The discretization uses the [`hydrostatic_reconstruction_chen_noelle`](@ref) on the conservative
variables.

Should be used together with [`Trixi.FluxHydrostaticReconstruction`](@extref) and
[`hydrostatic_reconstruction_chen_noelle`](@ref) in the surface flux to ensure consistency.

Further details on the hydrostatic reconstruction and its motivation can be found in
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI:10.1137/15M1053074](https://dx.doi.org/10.1137/15M1053074)
"""
@inline function flux_nonconservative_chen_noelle(u_ll, u_rr, orientation::Integer,
                                                  equations::ShallowWaterEquations2D)
    # Pull the water height and bottom topography on the left
    h_ll, _, _, b_ll = u_ll
    h_rr, _, _, b_rr = u_rr

    H_ll = h_ll + b_ll
    H_rr = h_rr + b_rr

    b_star = min(max(b_ll, b_rr), min(H_ll, H_rr))

    # Create the hydrostatic reconstruction for the left solution state
    u_ll_star, _ = hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, equations)

    # Copy the reconstructed water height for easier to read code
    h_ll_star = u_ll_star[1]

    z = zero(eltype(u_ll))

    g = equations.gravity
    if orientation == 1
        f = SVector(z,
                    -g * (h_ll_star + h_ll) * (b_ll - b_star),
                    z, z)
    else # orientation == 2
        f = SVector(z, z,
                    -g * (h_ll_star + h_ll) * (b_ll - b_star),
                    z)
    end

    return f
end

@inline function flux_nonconservative_chen_noelle(u_ll, u_rr,
                                                  normal_direction::AbstractVector,
                                                  equations::ShallowWaterEquations2D)
    # Pull the water height and bottom topography on the left
    h_ll, _, _, b_ll = u_ll
    h_rr, _, _, b_rr = u_rr

    H_ll = h_ll + b_ll
    H_rr = h_rr + b_rr

    b_star = min(max(b_ll, b_rr), min(H_ll, H_rr))

    # Create the hydrostatic reconstruction for the left solution state
    u_ll_star, _ = hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, equations)

    # Copy the reconstructed water height for easier to read code
    h_ll_star = u_ll_star[1]

    f2 = -normal_direction[1] * equations.gravity * (h_ll_star + h_ll) *
         (b_ll - b_star)
    f3 = -normal_direction[2] * equations.gravity * (h_ll_star + h_ll) *
         (b_ll - b_star)

    # First and last equations do not have a nonconservative flux
    f1 = f4 = zero(eltype(u_ll))

    return SVector(f1, f2, f3, f4)
end

"""
    flux_fjordholm_etal(u_ll, u_rr, orientation_or_normal_direction,
                        equations::ShallowWaterEquations2D)

Total energy conservative (mathematical entropy for shallow water equations). When the bottom topography
is nonzero this should only be used as a surface flux otherwise the scheme will not be well-balanced.
For well-balancedness in the volume flux use [`flux_wintermeyer_etal`](@ref).

Details are available in Eq. (4.1) in the paper:
- Ulrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011)
  Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography
  [DOI: 10.1016/j.jcp.2011.03.042](https://doi.org/10.1016/j.jcp.2011.03.042)
"""
@inline function Trixi.flux_fjordholm_etal(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquations2D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Average each factor of products in flux
    h_avg = 0.5f0 * (h_ll + h_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.25f0 * equations.gravity * (h_ll^2 + h_rr^2)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = h_avg * v1_avg
        f2 = f1 * v1_avg + p_avg
        f3 = f1 * v2_avg
    else
        f1 = h_avg * v2_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg + p_avg
    end

    return SVector(f1, f2, f3, 0)
end

@inline function Trixi.flux_fjordholm_etal(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::ShallowWaterEquations2D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Average each factor of products in flux
    h_avg = 0.5f0 * (h_ll + h_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    h2_avg = 0.5f0 * (h_ll^2 + h_rr^2)
    p_avg = 0.5f0 * equations.gravity * h2_avg
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)

    # Calculate fluxes depending on normal_direction
    f1 = h_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]

    return SVector(f1, f2, f3, 0)
end

"""
    flux_wintermeyer_etal(u_ll, u_rr, orientation_or_normal_direction,
                          equations::ShallowWaterEquations2D)

Total energy conservative (mathematical entropy for shallow water equations) split form.
When the bottom topography is nonzero this scheme will be well-balanced when used as a `volume_flux`.
For the `surface_flux` either [`flux_wintermeyer_etal`](@ref) or [`flux_fjordholm_etal`](@ref) can
be used to ensure well-balancedness and entropy conservation.

Further details are available in Theorem 1 of the paper:
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017)
  An entropy stable nodal discontinuous Galerkin method for the two dimensional
  shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry
  [DOI: 10.1016/j.jcp.2017.03.036](https://doi.org/10.1016/j.jcp.2017.03.036)
"""
@inline function Trixi.flux_wintermeyer_etal(u_ll, u_rr, orientation::Integer,
                                             equations::ShallowWaterEquations2D)
    # Unpack left and right state
    h_ll, h_v1_ll, h_v2_ll, _ = u_ll
    h_rr, h_v1_rr, h_v2_rr, _ = u_rr

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Average each factor of products in flux
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * equations.gravity * h_ll * h_rr

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = 0.5f0 * (h_v1_ll + h_v1_rr)
        f2 = f1 * v1_avg + p_avg
        f3 = f1 * v2_avg
    else
        f1 = 0.5f0 * (h_v2_ll + h_v2_rr)
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg + p_avg
    end

    return SVector(f1, f2, f3, 0)
end

@inline function Trixi.flux_wintermeyer_etal(u_ll, u_rr,
                                             normal_direction::AbstractVector,
                                             equations::ShallowWaterEquations2D)
    # Unpack left and right state
    h_ll, h_v1_ll, h_v2_ll, _ = u_ll
    h_rr, h_v1_rr, h_v2_rr, _ = u_rr

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Average each factor of products in flux
    h_v1_avg = 0.5f0 * (h_v1_ll + h_v1_rr)
    h_v2_avg = 0.5f0 * (h_v2_ll + h_v2_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * equations.gravity * h_ll * h_rr

    # Calculate fluxes depending on normal_direction
    f1 = h_v1_avg * normal_direction[1] + h_v2_avg * normal_direction[2]
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]

    return SVector(f1, f2, f3, 0)
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the bottom topography
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::ShallowWaterEquations2D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations)
    diss = -0.5f0 * λ * (u_rr - u_ll)
    return SVector(diss[1], diss[2], diss[3], 0)
end

# Specialized `FluxHLL` to avoid spurious dissipation in the bottom topography
@inline function (numflux::FluxHLL)(u_ll, u_rr, orientation_or_normal_direction,
                                    equations::ShallowWaterEquations2D)
    λ_min, λ_max = numflux.min_max_speed(u_ll, u_rr, orientation_or_normal_direction,
                                         equations)

    if λ_min >= 0 && λ_max >= 0
        return flux(u_ll, orientation_or_normal_direction, equations)
    elseif λ_max <= 0 && λ_min <= 0
        return flux(u_rr, orientation_or_normal_direction, equations)
    else
        f_ll = flux(u_ll, orientation_or_normal_direction, equations)
        f_rr = flux(u_rr, orientation_or_normal_direction, equations)
        inv_λ_max_minus_λ_min = inv(λ_max - λ_min)
        factor_ll = λ_max * inv_λ_max_minus_λ_min
        factor_rr = λ_min * inv_λ_max_minus_λ_min
        factor_diss = λ_min * λ_max * inv_λ_max_minus_λ_min
        diss = u_rr - u_ll
        return factor_ll * f_ll - factor_rr * f_rr +
               factor_diss * SVector(diss[1], diss[2], diss[3], zero(eltype(u_ll)))
    end
end

"""
    dissipation_roe(u_ll, u_rr, orientation_or_normal_direction,
                                    equations::ShallowWaterEquations2D)
Roe-type dissipation term for the [`ShallowWaterEquations2D`](@ref). To create the classical Roe solver,
this dissipation term can be combined with [`Trixi.flux_central`](@extref) using [`Trixi.FluxPlusDissipation`](@extref).

For details on the Roe linearization see Chapter 15.3.2 and Chapter 21.7 for the two-dimensional
shallow water equations of the book:
- Randall J. LeVeque (2002)
  Finite Volume Methods for Hyperbolic Problems
  [DOI: 10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253)
"""
@inline function dissipation_roe(u_ll, u_rr, normal_direction::AbstractVector,
                                 equations::ShallowWaterEquations2D)
    g = equations.gravity
    z = zero(eltype(u_ll))

    # Use the `normal_vector` to match how derived
    s_hat = Trixi.norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `s_hat` later
    normal = normal_direction / s_hat

    # Get velocities and waterheights
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Compute Roe averages
    h_avg = 0.5f0 * (h_ll + h_rr)
    v1_avg = (sqrt(h_ll) * v1_ll + sqrt(h_rr) * v1_rr) /
             (sqrt(h_ll) + sqrt(h_rr))
    v2_avg = (sqrt(h_ll) * v2_ll + sqrt(h_rr) * v2_rr) /
             (sqrt(h_ll) + sqrt(h_rr))
    c_avg = (sqrt(g * h_avg))
    vn_avg = normal[1] * v1_avg + normal[2] * v2_avg

    # Compute the eigenvalues
    λ1 = vn_avg - c_avg
    λ2 = vn_avg
    λ3 = vn_avg + c_avg

    # Eigenvector matrix
    r11 = 1
    r12 = 0
    r13 = 1

    r21 = v1_avg - c_avg * normal[1]
    r22 = -normal[2]
    r23 = v1_avg + c_avg * normal[1]

    r31 = v2_avg - c_avg * normal[2]
    r32 = normal[1]
    r33 = v2_avg + c_avg * normal[2]

    R = @SMatrix [[r11 r12 r13]; [r21 r22 r23]; [r31 r32 r33]]

    # Inverse eigenvector matrix
    inv_2c = inv(2 * c_avg)

    r11_inv = (vn_avg + c_avg) * inv_2c
    r12_inv = -normal[1] * inv_2c
    r13_inv = -normal[2] * inv_2c

    r21_inv = -(-normal[2] * v1_avg + normal[1] * v2_avg)
    r22_inv = -normal[2]
    r23_inv = normal[1]

    r31_inv = -(vn_avg - c_avg) * inv_2c
    r32_inv = normal[1] * inv_2c
    r33_inv = normal[2] * inv_2c

    R_inv = @SMatrix [[r11_inv r12_inv r13_inv]; [r21_inv r22_inv r23_inv];
                      [r31_inv r32_inv r33_inv]]

    # Eigenvalue absolute value matrix
    Λ_abs = @SMatrix [[abs(λ1) z z]; [z abs(λ2) z]; [z z abs(λ3)]]

    # Compute the jump in conserved variables, excluding the bottom topography
    u_jump = @views (u_rr - u_ll)[1:3]

    diss = SVector(-0.5f0 * R * Λ_abs * R_inv * u_jump)

    return SVector(diss[1], diss[2], diss[3], z) * s_hat
end

@inline function dissipation_roe(u_ll, u_rr, orientation::Integer,
                                 equations::ShallowWaterEquations2D)
    if orientation == 1
        return dissipation_roe(u_ll, u_rr, SVector(1, 0), equations)
    else # orientation == 2
        return dissipation_roe(u_ll, u_rr, SVector(0, 1), equations)
    end
end

"""
    min_max_speed_chen_noelle(u_ll, u_rr, orientation::Integer,
                              equations::ShallowWaterEquations2D)
    min_max_speed_chen_noelle(u_ll, u_rr, normal_direction::AbstractVector,
                              equations::ShallowWaterEquations2D)

Special estimate of the minimal and maximal wave speed of the shallow water equations for
the left and right states `u_ll, u_rr`. These approximate speeds are used for the HLL-type
numerical flux [`flux_hll_chen_noelle`](@ref). These wave speed estimates
together with a particular hydrostatic reconstruction technique guarantee
that the numerical flux is positive and satisfies an entropy inequality.

Further details on this hydrostatic reconstruction and its motivation can be found in
the reference below. The definition of the wave speeds are given in Equation (2.20).
- Guoxian Chen and Sebastian Noelle (2017)
  A new hydrostatic reconstruction scheme based on subcell reconstructions
  [DOI:10.1137/15M1053074](https://dx.doi.org/10.1137/15M1053074)
"""
@inline function min_max_speed_chen_noelle(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    a_ll = sqrt(equations.gravity * h_ll)
    a_rr = sqrt(equations.gravity * h_rr)

    if orientation == 1 # x-direction
        λ_min = min(v1_ll - a_ll, v1_rr - a_rr, zero(eltype(u_ll)))
        λ_max = max(v1_ll + a_ll, v1_rr + a_rr, zero(eltype(u_ll)))
    else # y-direction
        λ_min = min(v2_ll - a_ll, v2_rr - a_rr, zero(eltype(u_ll)))
        λ_max = max(v2_ll + a_ll, v2_rr + a_rr, zero(eltype(u_ll)))
    end

    return λ_min, λ_max
end

@inline function min_max_speed_chen_noelle(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    norm_ = norm(normal_direction)

    a_ll = sqrt(equations.gravity * h_ll) * norm_
    a_rr = sqrt(equations.gravity * h_rr) * norm_

    λ_min = min(v_normal_ll - a_ll, v_normal_rr - a_rr, zero(eltype(u_ll)))
    λ_max = max(v_normal_ll + a_ll, v_normal_rr + a_rr, zero(eltype(u_ll)))

    return λ_min, λ_max
end

# Wave speed estimates used in the CFL based time step selection
@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterEquations2D)
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)

    c = sqrt(equations.gravity * h)
    return abs(v1) + c, abs(v2) + c
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquations2D)
    # Get the velocity quantities in the appropriate direction
    if orientation == 1
        v_ll, _ = velocity(u_ll, equations)
        v_rr, _ = velocity(u_rr, equations)
    else
        _, v_ll = velocity(u_ll, equations)
        _, v_rr = velocity(u_rr, equations)
    end

    # Calculate the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    c_ll = sqrt(equations.gravity * h_ll)
    c_rr = sqrt(equations.gravity * h_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::ShallowWaterEquations2D)
    # Extract and compute the velocities in the normal direction
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    c_ll = sqrt(equations.gravity * h_ll)
    c_rr = sqrt(equations.gravity * h_rr)

    # The normal velocities are already scaled by the norm
    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function Trixi.max_abs_speed(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterEquations2D)
    # Get the velocity quantities in the appropriate direction
    if orientation == 1
        v_ll, _ = velocity(u_ll, equations)
        v_rr, _ = velocity(u_rr, equations)
    else
        _, v_ll = velocity(u_ll, equations)
        _, v_rr = velocity(u_rr, equations)
    end

    # Calculate the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    c_ll = sqrt(equations.gravity * h_ll)
    c_rr = sqrt(equations.gravity * h_rr)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function Trixi.max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::ShallowWaterEquations2D)
    # Extract and compute the velocities in the normal direction
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the wave celerity on the left and right
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    c_ll = sqrt(equations.gravity * h_ll)
    c_rr = sqrt(equations.gravity * h_rr)

    norm_ = norm(normal_direction)
    # The normal velocities are already scaled by the norm
    return max(abs(v_ll) + c_ll * norm_, abs(v_rr) + c_rr * norm_)
end

# Calculate estimates for minimum and maximum wave speeds for HLL-type fluxes
@inline function Trixi.min_max_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    if orientation == 1 # x-direction
        λ_min = v1_ll - sqrt(equations.gravity * h_ll)
        λ_max = v1_rr + sqrt(equations.gravity * h_rr)
    else # y-direction
        λ_min = v2_ll - sqrt(equations.gravity * h_ll)
        λ_max = v2_rr + sqrt(equations.gravity * h_rr)
    end

    return λ_min, λ_max
end

@inline function Trixi.min_max_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    norm_ = norm(normal_direction)
    # The v_normals are already scaled by the norm
    λ_min = v_normal_ll - sqrt(equations.gravity * h_ll) * norm_
    λ_max = v_normal_rr + sqrt(equations.gravity * h_rr) * norm_

    return λ_min, λ_max
end

# More refined estimates for minimum and maximum wave speeds for HLL-type fluxes
@inline function Trixi.min_max_speed_davis(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    c_ll = sqrt(equations.gravity * h_ll)
    c_rr = sqrt(equations.gravity * h_rr)

    if orientation == 1 # x-direction
        λ_min = min(v1_ll - c_ll, v1_rr - c_rr)
        λ_max = max(v1_ll + c_ll, v1_rr + c_rr)
    else # y-direction
        λ_min = min(v2_ll - c_ll, v2_rr - c_rr)
        λ_max = max(v2_ll + c_ll, v2_rr + c_rr)
    end

    return λ_min, λ_max
end

@inline function Trixi.min_max_speed_davis(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    norm_ = norm(normal_direction)
    c_ll = sqrt(equations.gravity * h_ll) * norm_
    c_rr = sqrt(equations.gravity * h_rr) * norm_

    v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # The v_normals are already scaled by the norm
    λ_min = min(v_normal_ll - c_ll, v_normal_rr - c_rr)
    λ_max = max(v_normal_ll + c_ll, v_normal_rr + c_rr)

    return λ_min, λ_max
end

@inline function Trixi.min_max_speed_einfeldt(u_ll, u_rr, orientation::Integer,
                                              equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    c_ll = sqrt(equations.gravity * h_ll)
    c_rr = sqrt(equations.gravity * h_rr)

    if orientation == 1 # x-direction
        v_roe, c_roe = calc_wavespeed_roe(u_ll, u_rr, orientation, equations)
        λ_min = min(v1_ll - c_ll, v_roe - c_roe)
        λ_max = max(v1_rr + c_rr, v_roe + c_roe)
    else # y-direction
        v_roe, c_roe = calc_wavespeed_roe(u_ll, u_rr, orientation, equations)
        λ_min = min(v2_ll - c_ll, v_roe - c_roe)
        λ_max = max(v2_rr + c_rr, v_roe + c_roe)
    end

    return λ_min, λ_max
end

@inline function Trixi.min_max_speed_einfeldt(u_ll, u_rr,
                                              normal_direction::AbstractVector,
                                              equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    norm_ = norm(normal_direction)

    c_ll = sqrt(equations.gravity * h_ll) * norm_
    c_rr = sqrt(equations.gravity * h_rr) * norm_

    v_normal_ll = (v1_ll * normal_direction[1] + v2_ll * normal_direction[2])
    v_normal_rr = (v1_rr * normal_direction[1] + v2_rr * normal_direction[2])

    v_roe, c_roe = calc_wavespeed_roe(u_ll, u_rr, normal_direction, equations)
    λ_min = min(v_normal_ll - c_ll, v_roe - c_roe)
    λ_max = max(v_normal_rr + c_rr, v_roe + c_roe)

    return λ_min, λ_max
end

"""
    calc_wavespeed_roe(u_ll, u_rr, direction::Integer,
                       equations::ShallowWaterEquations2D)

Calculate Roe-averaged velocity `v_roe` and wavespeed `c_roe = sqrt{g * h_roe}` depending on direction.
See for instance equation (62) in
- Paul A. Ullrich, Christiane Jablonowski, and Bram van Leer (2010)
  High-order finite-volume methods for the shallow-water equations on the sphere
  [DOI: 10.1016/j.jcp.2010.04.044](https://doi.org/10.1016/j.jcp.2010.04.044)
Or [this slides](https://faculty.washington.edu/rjl/classes/am574w2011/slides/am574lecture20nup3.pdf),
slides 8 and 9.
"""
@inline function calc_wavespeed_roe(u_ll, u_rr, orientation::Integer,
                                    equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    h_roe = 0.5f0 * (h_ll + h_rr)
    c_roe = sqrt(equations.gravity * h_roe)

    h_ll_sqrt = sqrt(h_ll)
    h_rr_sqrt = sqrt(h_rr)

    if orientation == 1 # x-direction
        v_roe = (h_ll_sqrt * v1_ll + h_rr_sqrt * v1_rr) / (h_ll_sqrt + h_rr_sqrt)
    else # y-direction
        v_roe = (h_ll_sqrt * v2_ll + h_rr_sqrt * v2_rr) / (h_ll_sqrt + h_rr_sqrt)
    end

    return v_roe, c_roe
end

@inline function calc_wavespeed_roe(u_ll, u_rr, normal_direction::AbstractVector,
                                    equations::ShallowWaterEquations2D)
    h_ll = waterheight(u_ll, equations)
    v1_ll, v2_ll = velocity(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    norm_ = norm(normal_direction)

    h_roe = 0.5f0 * (h_ll + h_rr)
    c_roe = sqrt(equations.gravity * h_roe) * norm_

    h_ll_sqrt = sqrt(h_ll)
    h_rr_sqrt = sqrt(h_rr)

    v1_roe = (h_ll_sqrt * v1_ll + h_rr_sqrt * v1_rr) / (h_ll_sqrt + h_rr_sqrt)
    v2_roe = (h_ll_sqrt * v2_ll + h_rr_sqrt * v2_rr) / (h_ll_sqrt + h_rr_sqrt)

    v_roe = (v1_roe * normal_direction[1] + v2_roe * normal_direction[2])

    return v_roe, c_roe
end

@inline function Trixi.rotate_to_x(u, normal_vector,
                                   equations::ShallowWaterEquations2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    # Apply the 2D rotation matrix with normal and tangent directions of the form
    # [ 1    0    0   0;
    #   0   n_1  n_2  0;
    #   0   t_1  t_2  0;
    #   0    0    0   1 ]
    # where t_1 = -n_2 and t_2 = n_1

    return SVector(u[1],
                   c * u[2] + s * u[3],
                   -s * u[2] + c * u[3],
                   u[4])
end

@inline function Trixi.rotate_from_x(u, normal_vector,
                                     equations::ShallowWaterEquations2D)
    # cos and sin of the angle between the x-axis and the normalized normal_vector are
    # the normalized vector's x and y coordinates respectively (see unit circle).
    c = normal_vector[1]
    s = normal_vector[2]

    # Apply the 2D back-rotation matrix with normal and tangent directions of the form
    # [ 1    0    0   0;
    #   0   n_1  t_1  0;
    #   0   n_2  t_2  0;
    #   0    0    0   1 ]
    # where t_1 = -n_2 and t_2 = n_1

    return SVector(u[1],
                   c * u[2] - s * u[3],
                   s * u[2] + c * u[3],
                   u[4])
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterEquations2D)
    h, h_v1, h_v2, _ = u

    v1 = h_v1 / h
    v2 = h_v2 / h
    return SVector(v1, v2)
end

@inline function Trixi.velocity(u, orientation::Int, equations::ShallowWaterEquations2D)
    h = u[1]
    v = u[orientation + 1] / h
    return v
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterEquations2D)
    h, _, _, b = u

    H = h + b
    v1, v2 = velocity(u, equations)
    return SVector(H, v1, v2, b)
end

# Convert conservative variables to entropy
# Note, only the first three are the entropy variables, the fourth entry still
# just carries the bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::ShallowWaterEquations2D)
    h, _, _, b = u

    v1, v2 = velocity(u, equations)
    v_square = v1^2 + v2^2

    w1 = equations.gravity * (h + b) - 0.5f0 * v_square
    w2 = v1
    w3 = v2
    return SVector(w1, w2, w3, b)
end

# Convert entropy variables to conservative
@inline function Trixi.entropy2cons(w, equations::ShallowWaterEquations2D)
    w1, w2, w3, b = w

    h = (w1 + 0.5f0 * (w2^2 + w3^2)) / equations.gravity - b
    h_v1 = h * w2
    h_v2 = h * w3
    return SVector(h, h_v1, h_v2, b)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterEquations2D)
    H, v1, v2, b = prim

    h = H - b
    h_v1 = h * v1
    h_v2 = h * v2
    return SVector(h, h_v1, h_v2, b)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterEquations2D)
    return u[1]
end

@inline function Trixi.pressure(u, equations::ShallowWaterEquations2D)
    h = waterheight(u, equations)
    p = 0.5f0 * equations.gravity * h^2
    return p
end

@inline function Trixi.waterheight_pressure(u, equations::ShallowWaterEquations2D)
    return waterheight(u, equations) * Trixi.pressure(u, equations)
end

# Entropy function for the shallow water equations is the total energy
@inline function Trixi.entropy(cons, equations::ShallowWaterEquations2D)
    Trixi.energy_total(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline function Trixi.energy_total(cons, equations::ShallowWaterEquations2D)
    h, h_v1, h_v2, b = cons

    e = (h_v1^2 + h_v2^2) / (2 * h) + 0.5f0 * equations.gravity * h^2 +
        equations.gravity * h * b
    return e
end

# Calculate kinetic energy for a conservative state `cons`
@inline function Trixi.energy_kinetic(u, equations::ShallowWaterEquations2D)
    h, h_v1, h_v2, _ = u
    return (h_v1^2 + h_v2^2) / (2 * h)
end

# Calculate potential energy for a conservative state `cons`
@inline function Trixi.energy_internal(cons, equations::ShallowWaterEquations2D)
    return energy_total(cons, equations) - energy_kinetic(cons, equations)
end

# Calculate the error for the "lake-at-rest" test case where H = h+b should
# be a constant value over time. Note, assumes there is a single reference
# water height `H0` with which to compare.
@inline function Trixi.lake_at_rest_error(u, equations::ShallowWaterEquations2D)
    h, _, _, b = u

    # For well-balancedness testing with possible wet/dry regions the reference
    # water height `H0` accounts for the possibility that the bottom topography
    # can emerge out of the water as well as for the threshold offset to avoid
    # division by a "hard" zero water heights as well.
    H0_wet_dry = max(equations.H0, b + equations.threshold_limiter)

    return abs(H0_wet_dry - (h + b))
end
end # @muladd
