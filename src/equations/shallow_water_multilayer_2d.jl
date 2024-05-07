# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    ShallowWaterMultiLayerEquations2D(gravity, H0, rhos)

Multi-Layer Shallow Water equations (MLSWE) in two space dimension. The equations are given by
```math
\left\{
	\begin{aligned}			
		&\partial_t h_m + \partial_x h_mv_m = 0,\\
		&\partial h_mv1_m + \partial_x h_mv1_m^2 + \partial_y h_mv1_mv2_m = -gh_m\partial_x \bigg(b + \sum\limits_{k\geq j}h_k + \sum\limits_{k<m}\frac{\rho_k}{\rho_m}h_k \bigg)
        &\partial h_mv2_m + \partial_x h_mv1_mv2_m + \partial_y h_mv2_m^2 = -gh_m\partial_y \bigg(b + \sum\limits_{k\geq j}h_k + \sum\limits_{k<m}\frac{\rho_k}{\rho_m}h_k \bigg)
	\end{aligned}
\right.
```

where ``m = 1, 2, ..., M`` is the layer index and the unknown variables are the water height ``h`` and
the velocities ``v1, v2`` in both spatial dimensions .  Furthermore, ``g`` denotes the gravitational 
constant, ``b(x)`` the bottom topography and ``\rho_m`` the m-th layer density, that must be chosen such that 
``\rho_1 < \rho_2 < ... < \rho_M``, to ensure that different layers are ordered from top to bottom, with 
increasing density.

We use a specific formulation of the system, where the pressure term is reformulated as a 
nonconservative term, which has some benefits for the design of well-balanced schemes.

The additional quantity ``H_0`` is also available to store a reference value for the total water
height that is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

Also, there are two thresholds which prevent numerical problems as well as instabilities. The limiters are 
used in [`PositivityPreservingLimiterShallowWater`](@ref) on the water height. `threshold_limiter` 
acts as a (small) shift on the initial condition and cutoff before the next time step, whereas 
`threshold_desingularization` is used in the velocity desingularization. A third 
`threshold_partially_wet` is applied on the water height to define "partially wet" elements in 
[`IndicatorHennemannGassnerShallowWater`](@ref), that are then calculated with a pure FV method to
ensure well-balancedness. For `Float64` no threshold needs to be passed, as default values are 
defined within the struct. For other number formats `threshold_desingularization` and `threshold_partially_wet` 
must be provided.

The bottom topography function ``b(x)`` is set inside the initial condition routine
for a particular problem setup.

In addition to the unknowns, Trixi currently stores the bottom topography values at the
approximation points despite being fixed in time. This is done for convenience of computing the
bottom topography gradients on the fly during the approximation as well as computing auxiliary
quantities like the total water height ``H`` or the entropy variables.
This affects the implementation and use of these equations in various ways:
* The flux values corresponding to the bottom topography must be zero.
* The bottom topography values must be included when defining initial conditions, boundary
  conditions or source terms.
* [`AnalysisCallback`](@ref) analyzes this variable.
* Trixi's visualization tools will visualize the bottom topography by default.

A good introduction for the MLSWE is available in Chapter 12 of the book:
    - Benoit Cushman-Roisin (2011)\
      Introduction to geophyiscal fluid dynamics: physical and numerical aspects\
      <https://www.sciencedirect.com/bookseries/international-geophysics/vol/101/suppl/C>\
      ISBN: 978-0-12-088759-0
"""
struct ShallowWaterMultiLayerEquations2D{NVARS, NLAYERS, RealT <: Real} <:
       AbstractShallowWaterMultiLayerEquations{2, NVARS, NLAYERS}
    gravity::RealT   # gravitational constant
    H0::RealT        # constant "lake-at-rest" total water height
    threshold_limiter::RealT    # threshold for the positivity-limiter
    threshold_desingularization::RealT  # threshold for velocity desingularization
    threshold_partially_wet::RealT  # threshold to define partially wet elements
    rhos::SVector{NLAYERS, RealT} # Vector of layer densities

    function ShallowWaterMultiLayerEquations2D{NVARS, NLAYERS, RealT}(gravity::RealT,
                                                                      H0::RealT,
                                                                      threshold_limiter::RealT,
                                                                      threshold_desingularization::RealT,
                                                                      threshold_partially_wet::RealT,
                                                                      rhos::SVector{NLAYERS,
                                                                                    RealT}) where {
                                                                                                   NVARS,
                                                                                                   NLAYERS,
                                                                                                   RealT <:
                                                                                                   Real
                                                                                                   }
        # Ensure that layer densities are all positive and in increasing order
        issorted(rhos) ||
            throw(ArgumentError("densities must be in increasing order (rhos[1] < rhos[2] < ... < rhos[NLAYERS])"))
        min(rhos...) > 0 || throw(ArgumentError("densities must be positive"))

        new(gravity, H0, threshold_limiter, threshold_desingularization,
            threshold_partially_wet, rhos)
    end
end

# Allow for flexibility to set the gravitational constant within an elixir depending on the
# application where `gravity_constant=1.0` or `gravity_constant=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases. 
function ShallowWaterMultiLayerEquations2D(; gravity_constant,
                                           H0 = zero(gravity_constant),
                                           threshold_limiter = nothing,
                                           threshold_desingularization = nothing,
                                           threshold_partially_wet = nothing,
                                           rhos)

    # Promote all variables to a common type
    _rhos = promote(rhos...)
    RealT = promote_type(eltype(_rhos), eltype(gravity_constant), eltype(H0))
    __rhos = SVector(map(RealT, _rhos))

    # Set default values for thresholds
    if threshold_limiter === nothing
        threshold_limiter = 5 * eps(RealT)
    end
    if threshold_desingularization === nothing
        threshold_desingularization = default_threshold_desingularization(RealT)
    end
    if threshold_partially_wet === nothing
        threshold_partially_wet = default_threshold_partially_wet(RealT)
    end

    # Extract number of layers and variables
    NLAYERS = length(rhos)
    NVARS = 3 * NLAYERS + 1

    return ShallowWaterMultiLayerEquations2D{NVARS, NLAYERS, RealT}(gravity_constant,
                                                                    H0,
                                                                    threshold_limiter,
                                                                    threshold_desingularization,
                                                                    threshold_partially_wet,
                                                                    __rhos)
end

@inline function Base.real(::ShallowWaterMultiLayerEquations2D{NVARS, NLAYERS, RealT}) where {
                                                                                              NVARS,
                                                                                              NLAYERS,
                                                                                              RealT <:
                                                                                              Real
                                                                                              }
    RealT
end

Trixi.have_nonconservative_terms(::ShallowWaterMultiLayerEquations2D) = True()

function Trixi.varnames(::typeof(cons2cons),
                        equations::ShallowWaterMultiLayerEquations2D)
    waterheight = ntuple(n -> "h" * string(n), Val(nlayers(equations)))
    momentum_1 = ntuple(n -> "h" * string(n) * "_v1", Val(nlayers(equations)))
    momentum_2 = ntuple(n -> "h" * string(n) * "_v2", Val(nlayers(equations)))
    return (waterheight..., momentum_1..., momentum_2..., "b")
end

# We use the total layer heights, H = ∑h + b as primitive variables for easier visualization and setting initial
# conditions
function Trixi.varnames(::typeof(cons2prim),
                        equations::ShallowWaterMultiLayerEquations2D)
    total_layer_height = ntuple(n -> "H" * string(n), Val(nlayers(equations)))
    velocity_1 = ntuple(n -> "v" * string(n) * "_1", Val(nlayers(equations)))
    velocity_2 = ntuple(n -> "v" * string(n) * "_2", Val(nlayers(equations)))
    return (total_layer_height..., velocity_1..., velocity_2..., "b")
end

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterMultiLayerEquations2D)

A smooth initial condition for a three-layer configuration used for convergence tests in combination with
[`source_terms_convergence_test`](@ref) (and
[`Trixi.BoundaryConditionDirichlet`](@extref) in non-periodic domains).
"""
function Trixi.initial_condition_convergence_test(x, t,
                                                  equations::ShallowWaterMultiLayerEquations2D)
    # Check if the equation has been set up with three layers
    if nlayers(equations) != 3
        throw(ArgumentError("This initial condition is only valid for three layers"))
    end

    # Some constants are chosen such that the function is periodic on the domain [0,sqrt(2)]
    ω = pi * sqrt(2.0)

    H3 = 1.5 + 0.1 * cos(ω * x[1] + t) + 0.1 * cos(ω * x[2] + t)
    H2 = 2.0 + 0.1 * sin(ω * x[1] + t) + 0.1 * sin(ω * x[2] + t)
    H1 = 4.0 + 0.1 * cos(ω * x[1] + t) + 0.1 * cos(ω * x[2] + t)
    H = (H1, H2, H3)
    v1 = (0.8, 0.8, 0.8)
    v2 = (1.0, 1.0, 1.0)

    b = 1.0 + 0.1 * cos(ω * x[1]) + 0.1 * cos(ω * x[2])

    return prim2cons(SVector(H..., v1..., v2..., b), equations)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterMultiLayerEquations2D)

Source terms used for convergence tests with a three-layer configuration in combination with
[`initial_condition_convergence_test`](@ref)
(and [`Trixi.BoundaryConditionDirichlet`](@extref)
in non-periodic domains).
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterMultiLayerEquations2D)
    # Same settings as in `initial_condition_convergence_test`. Some derivative simplify because
    # this manufactured solution velocity is taken to be constant
    ω = pi * sqrt(2.0)
    g = equations.gravity

    du1 = (-0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
           0.1 * sin(t + x[2] * ω) -
           0.1 * cos(t + x[1] * ω) - 0.1 * cos(t + x[2] * ω) * ω +
           0.8 * (-0.1 * sin(t + x[1] * ω) * ω - 0.1 * cos(t + x[1] * ω) * ω) -
           0.1 * sin(t + x[2] * ω) * ω)

    du2 = (0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) + 0.1 * sin(t + x[2] * ω) +
           0.1 * cos(t + x[1] * ω) + 0.1 * cos(t + x[2] * ω) * ω +
           0.8 * (0.1 * sin(t + x[1] * ω) * ω + 0.1 * cos(t + x[1] * ω) * ω) +
           0.1 * sin(t + x[2] * ω) * ω)

    du3 = (-0.1 * sin(t + x[1] * ω) - 0.1 * sin(t + x[2] * ω) +
           0.1 * sin(x[2] * ω) * ω +
           0.8 * (0.1 * sin(x[1] * ω) * ω - 0.1 * sin(t + x[1] * ω) * ω) -
           0.1 * sin(t + x[2] * ω) * ω)

    du4 = (0.8 * (-0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
            0.1 * sin(t + x[2] * ω) -
            0.1 * cos(t + x[1] * ω)) +
           0.8 * (-0.1 * cos(t + x[2] * ω) * ω - 0.1 * sin(t + x[2] * ω) * ω) +
           0.8^2 * (-0.1 * sin(t + x[1] * ω) * ω - 0.1 * cos(t + x[1] * ω) * ω) +
           g *
           (2.0 + 0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
            0.1 * sin(t + x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) *
           (-0.1 * sin(t + x[1] * ω) * ω - 0.1 * cos(t + x[1] * ω) * ω) +
           0.1 * g *
           (2.0 + 0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
            0.1 * sin(t + x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) * cos(t + x[1] * ω) * ω)

    du5 = (0.8 * (0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) +
            0.1 * sin(t + x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) +
           0.8 * (0.1 * cos(t + x[2] * ω) * ω + 0.1 * sin(t + x[2] * ω) * ω) +
           0.8^2 * (0.1 * sin(t + x[1] * ω) * ω + 0.1 * cos(t + x[1] * ω) * ω) +
           g *
           (0.5 - 0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) +
            0.1 * sin(t + x[2] * ω) -
            0.1 * cos(t + x[1] * ω)) *
           (-0.1 * sin(t + x[1] * ω) * ω +
            0.9 * (-0.1 * sin(t + x[1] * ω) * ω - 0.1 * cos(t + x[1] * ω) * ω)) +
           g *
           (0.5 - 0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) +
            0.1 * sin(t + x[2] * ω) -
            0.1 * cos(t + x[1] * ω)) *
           (0.1 * sin(t + x[1] * ω) * ω + 0.1 * cos(t + x[1] * ω) * ω))

    du6 = (0.8 * (-0.1 * sin(t + x[1] * ω) - 0.1 * sin(t + x[2] * ω)) +
           0.8 * (0.1 * sin(x[2] * ω) * ω - 0.1 * sin(t + x[2] * ω) * ω) +
           0.8^2 * (0.1 * sin(x[1] * ω) * ω - 0.1 * sin(t + x[1] * ω) * ω) +
           g *
           (0.5 - 0.1 * cos(x[1] * ω) + 0.1 * cos(t + x[2] * ω) - 0.1 * cos(x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) *
           (0.1 * sin(x[1] * ω) * ω - 0.1 * sin(t + x[1] * ω) * ω) +
           g *
           (0.5 - 0.1 * cos(x[1] * ω) + 0.1 * cos(t + x[2] * ω) - 0.1 * cos(x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) *
           (-0.1 * sin(x[1] * ω) * ω +
            0.9 / 1.1 * (-0.1 * sin(t + x[1] * ω) * ω - 0.1 * cos(t + x[1] * ω) * ω) +
            1.0 / 1.1 * (0.1 * sin(t + x[1] * ω) * ω + 0.1 * cos(t + x[1] * ω) * ω)))

    du7 = (-0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
           0.1 * sin(t + x[2] * ω) -
           0.1 * cos(t + x[1] * ω) - 0.1 * cos(t + x[2] * ω) * ω +
           0.8 * (-0.1 * sin(t + x[1] * ω) * ω - 0.1 * cos(t + x[1] * ω) * ω) -
           0.1 * sin(t + x[2] * ω) * ω +
           0.1 * g *
           (2.0 + 0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
            0.1 * sin(t + x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) * cos(t + x[2] * ω) * ω +
           g *
           (2.0 + 0.1 * cos(t + x[2] * ω) - 0.1 * sin(t + x[1] * ω) -
            0.1 * sin(t + x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) *
           (-0.1 * cos(t + x[2] * ω) * ω - 0.1 * sin(t + x[2] * ω) * ω))

    du8 = (0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) + 0.1 * sin(t + x[2] * ω) +
           0.1 * cos(t + x[1] * ω) + 0.1 * cos(t + x[2] * ω) * ω +
           0.8 * (0.1 * sin(t + x[1] * ω) * ω + 0.1 * cos(t + x[1] * ω) * ω) +
           0.1 * sin(t + x[2] * ω) * ω +
           g *
           (0.5 - 0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) +
            0.1 * sin(t + x[2] * ω) -
            0.1 * cos(t + x[1] * ω)) *
           (0.1 * cos(t + x[2] * ω) * ω + 0.1 * sin(t + x[2] * ω) * ω) +
           g *
           (0.5 - 0.1 * cos(t + x[2] * ω) + 0.1 * sin(t + x[1] * ω) +
            0.1 * sin(t + x[2] * ω) -
            0.1 * cos(t + x[1] * ω)) *
           (0.9 * (-0.1 * cos(t + x[2] * ω) * ω - 0.1 * sin(t + x[2] * ω) * ω) -
            0.1 * sin(t + x[2] * ω) * ω))

    du9 = (-0.1 * sin(t + x[1] * ω) - 0.1 * sin(t + x[2] * ω) +
           0.1 * sin(x[2] * ω) * ω +
           0.8 * (0.1 * sin(x[1] * ω) * ω - 0.1 * sin(t + x[1] * ω) * ω) -
           0.1 * sin(t + x[2] * ω) * ω +
           g *
           (0.5 - 0.1 * cos(x[1] * ω) + 0.1 * cos(t + x[2] * ω) - 0.1 * cos(x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) *
           (0.9 / 1.1 * (-0.1 * cos(t + x[2] * ω) * ω - 0.1 * sin(t + x[2] * ω) * ω) +
            1.0 / 1.1 * (0.1 * cos(t + x[2] * ω) * ω + 0.1 * sin(t + x[2] * ω) * ω) -
            0.1 * sin(x[2] * ω) * ω) +
           g *
           (0.5 - 0.1 * cos(x[1] * ω) + 0.1 * cos(t + x[2] * ω) - 0.1 * cos(x[2] * ω) +
            0.1 * cos(t + x[1] * ω)) *
           (0.1 * sin(x[2] * ω) * ω - 0.1 * sin(t + x[2] * ω) * ω))

    return SVector(du1, du2, du3, du4, du5, du6, du7, du8, du9, zero(eltype(u)))
end

"""
    boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,
                                 equations::ShallowWaterMultiLayerEquations2D)
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
                                                    surface_flux_function,
                                                    equations::ShallowWaterMultiLayerEquations2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # Extract internal values
    h = waterheight(u_inner, equations)
    h_v1, h_v2 = momentum(u_inner, equations)
    b = u_inner[end]

    # compute the normal velocity
    u_normal = normal[1] * h_v1 + normal[2] * h_v2

    # create the "external" boundary solution state
    u_boundary = SVector(h...,
                         (h_v1 - 2.0 * u_normal * normal[1])...,
                         (h_v2 - 2.0 * u_normal * normal[2])...,
                         b)

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

    return flux
end

"""
    boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                 surface_flux_function, equations::ShallowWaterMultiLayerEquations2D)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation,
                                                    direction, x, t,
                                                    surface_flux_function,
                                                    equations::ShallowWaterMultiLayerEquations2D)
    # Extract internal values
    h = waterheight(u_inner, equations)
    h_v1, h_v2 = momentum(u_inner, equations)
    b = u_inner[end]

    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(h..., -h_v1..., h_v2..., b)
    else # orientation == 2
        u_boundary = SVector(h..., h_v1..., -h_v2..., b)
    end

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

# Calculate 2D advective portion of the flux for a single point
# Note, the bottom topography has no flux
@inline function Trixi.flux(u, orientation::Integer,
                            equations::ShallowWaterMultiLayerEquations2D)

    # Extract waterheight and momentum and compute velocities
    h_v1, h_v2 = momentum(u, equations)
    v1, v2 = velocity(u, equations)

    # Initialize flux vector
    f = zero(MVector{3 * nlayers(equations) + 1, real(equations)})

    # Calculate fluxes in each layer
    # The momentum flux simplifies as the pressure is included in the nonconservative term
    if orientation == 1
        for i in eachlayer(equations)
            f_h = h_v1[i]
            f_hv1 = f_h * v1[i]
            f_hv2 = f_h * v2[i]
            setlayer!(f, f_h, f_hv1, f_hv2, i, equations)
        end
    else # orientation == 2
        for i in eachlayer(equations)
            f_h = h_v2[i]
            f_hv1 = f_h * v1[i]
            f_hv2 = f_h * v2[i]
            setlayer!(f, f_h, f_hv1, f_hv2, i, equations)
        end
    end

    return SVector(f)
end

# Calculate 2D flux for a single point in the normal direction
# Note, this directional vector is not normalized and the bottom topography has no flux
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::ShallowWaterMultiLayerEquations2D)
    # Extract waterheight and momentum and compute velocities
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    h_v_normal = h .* v_normal

    # Initialize flux vector
    f = zero(MVector{3 * nlayers(equations) + 1, real(equations)})

    # Calculate fluxes in each layer
    # The momentum flux simplifies as the pressure is included in the nonconservative term
    for i in eachlayer(equations)
        f_h = h_v_normal[i]
        f_hv1 = f_h * v1[i]
        f_hv2 = f_h * v2[i]
        setlayer!(f, f_h, f_hv1, f_hv2, i, equations)
    end

    return SVector(f)
end

"""
    flux_nonconservative_ersing_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterMultiLayerEquations2D)
    flux_nonconservative_ersing_etal(u_ll, u_rr,
                                     normal_direction_ll::AbstractVector,
                                     normal_direction_average::AbstractVector,
                                     equations::ShallowWaterMultiLayerEquations2D)

Non-symmetric path-conservative two-point flux discretizing the nonconservative (source) term
that contains the gradients of the bottom topography and waterheights from the coupling between layers
and the nonconservative pressure formulation [`ShallowWaterMultiLayerEquations2D`](@ref).

When the bottom topography is nonzero this scheme will be well-balanced when used with [`flux_ersing_etal`](@ref).

In the two-layer setting this combination is equivalent to the fluxes in:
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on 
  curvilinear meshes
  [DOI: 10.1007/s10915-024-02451-2](https://doi.org/10.1007/s10915-024-02451-2)
"""
@inline function Trixi.flux_nonconservative_ersing_etal(u_ll, u_rr,
                                                        orientation::Integer,
                                                        equations::ShallowWaterMultiLayerEquations2D)
    # Pull the necessary left and right state information
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    b_rr = u_rr[end]
    b_ll = u_ll[end]

    # Compute the jumps
    h_jump = h_rr - h_ll
    b_jump = b_rr - b_ll
    g = equations.gravity

    # Initialize flux vector
    f = zero(MVector{3 * nlayers(equations) + 1, real(equations)})

    # Compute the nonconservative flux in each layer
    # where f_hv[i] = g * h[i] * (b + ∑h[k] + ∑σ[k] * h[k])_x and σ[k] = ρ[k] / ρ[i] denotes the 
    # density ratio of different layers
    for i in eachlayer(equations)
        f_hv = g * h_ll[i] * b_jump
        for j in eachlayer(equations)
            if j < i
                f_hv += g * h_ll[i] *
                        (equations.rhos[j] / equations.rhos[i] * h_jump[j])
            else # (i<j<nlayers) nonconservative formulation of the pressure
                f_hv += g * h_ll[i] * h_jump[j]
            end
        end

        if orientation == 1
            setindex!(f, f_hv, i + nlayers(equations))
        else # orientation == 2
            setindex!(f, f_hv, i + 2 * nlayers(equations))
        end
    end

    return SVector(f)
end

@inline function Trixi.flux_nonconservative_ersing_etal(u_ll, u_rr,
                                                        normal_direction_ll::AbstractVector,
                                                        normal_direction_average::AbstractVector,
                                                        equations::ShallowWaterMultiLayerEquations2D)
    # Pull the necessary left and right state information
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    b_rr = u_rr[end]
    b_ll = u_ll[end]

    # Compute the jumps
    h_jump = h_rr - h_ll
    b_jump = b_rr - b_ll
    g = equations.gravity

    # Initialize flux vector
    f = zero(MVector{3 * nlayers(equations) + 1, real(equations)})

    for i in eachlayer(equations)
        f_h = zero(real(equations))
        f_hv = g * h_ll[i] * b_jump
        for j in eachlayer(equations)
            if j < i
                f_hv += g * h_ll[i] *
                        (equations.rhos[j] / equations.rhos[i] * h_jump[j])
            else # (i<j<nlayers) nonconservative formulation of the pressure
                f_hv += g * h_ll[i] * h_jump[j]
            end
        end
        setlayer!(f, f_h, f_hv * normal_direction_average[1],
                  f_hv * normal_direction_average[2], i, equations)
    end

    return SVector(f)
end

"""
    flux_ersing_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterMultiLayerEquations2D)
    flux_ersing_etal(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::ShallowWaterMultiLayerEquations2D)

Total energy conservative (mathematical entropy for MLSWE) split form,
without the hydrostatic pressure.
When the bottom topography is nonzero this scheme will be well-balanced when used with the 
nonconservative [`flux_nonconservative_ersing_etal`](@ref).

To obtain an entropy stable formulation the `surface_flux` can be set as
`FluxPlusDissipation(flux_ersing_etal, DissipationLocalLaxFriedrichs()), flux_nonconservative_ersing_etal`.

In the two-layer setting this combination is equivalent to the fluxes in:
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on 
  curvilinear meshes
  [DOI: 10.1007/s10915-024-02451-2](https://doi.org/10.1007/s10915-024-02451-2)
"""
@inline function flux_ersing_etal(u_ll, u_rr,
                                  orientation::Integer,
                                  equations::ShallowWaterMultiLayerEquations2D)
    # Unpack left and right state
    h_v1_ll, h_v2_ll = momentum(u_ll, equations)
    h_v1_rr, h_v2_rr = momentum(u_rr, equations)

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Initialize flux vector
    f = zero(MVector{3 * nlayers(equations) + 1, real(equations)})

    # Calculate fluxes in each layer
    for i in eachlayer(equations)
        # Compute averages
        v1_avg = 0.5 * (v1_ll[i] + v1_rr[i])
        v2_avg = 0.5 * (v2_ll[i] + v2_rr[i])
        h_v1_avg = 0.5 * (h_v1_ll[i] + h_v1_rr[i])
        h_v2_avg = 0.5 * (h_v2_ll[i] + h_v2_rr[i])

        # Compute fluxes
        # The momentum flux simplifies as the pressure is included in the nonconservative term.
        if orientation == 1
            f_h = h_v1_avg
            f_hv1 = f_h * v1_avg
            f_hv2 = f_h * v2_avg
        else # orientation == 2
            f_h = h_v2_avg
            f_hv1 = f_h * v1_avg
            f_hv2 = f_h * v2_avg
        end
        setlayer!(f, f_h, f_hv1, f_hv2, i, equations)
    end

    return SVector(f)
end

@inline function flux_ersing_etal(u_ll, u_rr,
                                  normal_direction::AbstractVector,
                                  equations::ShallowWaterMultiLayerEquations2D)
    # Unpack left and right state
    h_v1_ll, h_v2_ll = momentum(u_ll, equations)
    h_v1_rr, h_v2_rr = momentum(u_rr, equations)

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Initialize flux vector
    f = zero(MVector{3 * nlayers(equations) + 1, real(equations)})

    # Calculate fluxes in each layer
    for i in eachlayer(equations)
        # Compute averages
        v1_avg = 0.5 * (v1_ll[i] + v1_rr[i])
        v2_avg = 0.5 * (v2_ll[i] + v2_rr[i])
        h_v1_avg = 0.5 * (h_v1_ll[i] + h_v1_rr[i])
        h_v2_avg = 0.5 * (h_v2_ll[i] + h_v2_rr[i])

        # Compute fluxes
        # The momentum flux simplifies as the pressure is included in the nonconservative term.
        f_h = h_v1_avg * normal_direction[1] + h_v2_avg * normal_direction[2]
        f_hv1 = f_h * v1_avg
        f_hv2 = f_h * v2_avg

        setlayer!(f, f_h, f_hv1, f_hv2, i, equations)
    end

    return SVector(f)
end

"""
    hydrostatic_reconstruction_ersing_etal(u_ll, u_rr, equations::ShallowWaterMultiLayerEquations2D)

A particular type of hydrostatic reconstruction of the water height and bottom topography to 
guarantee well-balancedness in the presence of wet/dry transitions and entropy stability for the 
[`ShallowWaterMultiLayerEquations2D`](@ref). 
The reconstructed solution states `u_ll_star` and `u_rr_star` variables are used to evaluate the 
surface numerical flux at the interface. The key idea is a piecewise linear reconstruction of the 
bottom topography and water height interfaces using subcells, where the bottom topography is allowed 
to be discontinuous. 
Use in combination with the generic numerical flux routine [`Trixi.FluxHydrostaticReconstruction`](@extref).

!!! warning "Experimental code"
    This is an experimental feature and may change in future releases.
"""
@inline function hydrostatic_reconstruction_ersing_etal(u_ll, u_rr,
                                                        equations::ShallowWaterMultiLayerEquations2D)
    # Unpack waterheight and bottom topographies
    h_ll = MVector(waterheight(u_ll, equations))
    h_rr = MVector(waterheight(u_rr, equations))
    b_ll = u_ll[end]
    b_rr = u_rr[end]

    # Get the velocities on either side
    v1_ll = MVector(velocity(u_ll, equations)[1])
    v2_ll = MVector(velocity(u_ll, equations)[2])
    v1_rr = MVector(velocity(u_rr, equations)[1])
    v2_rr = MVector(velocity(u_rr, equations)[2])

    threshold = equations.threshold_limiter

    # Ensure zero velocity at dry states
    for i in eachlayer(equations)
        if h_ll[i] <= threshold
            v1_ll[i] = zero(eltype(u_ll))
            v2_ll[i] = zero(eltype(u_ll))
        end
        if h_rr[i] <= threshold
            v1_rr[i] = zero(eltype(u_rr))
            v2_rr[i] = zero(eltype(u_rr))
        end
    end

    # Calculate total layer heights
    H_ll = waterheight(cons2prim(u_ll, equations), equations)
    H_rr = waterheight(cons2prim(u_rr, equations), equations)

    # Discontinuous reconstruction of the bottom topography
    b_ll_star = min(H_ll[1], max(b_rr, b_ll))
    b_rr_star = min(H_rr[1], max(b_rr, b_ll))

    # Calculate reconstructed total layer heights
    H_ll_star = max.(H_ll, b_ll_star)
    H_rr_star = max.(H_rr, b_rr_star)

    # Initialize reconstructed waterheights
    h_ll_star = zero(MVector{nlayers(equations), real(equations)})
    h_rr_star = zero(MVector{nlayers(equations), real(equations)})

    # Reconstruct the waterheights
    for i in eachlayer(equations)
        if i == nlayers(equations) # The lowest layer is measured from the bottom topography
            h_ll_star[i] = max(H_ll_star[i] - b_ll_star, threshold)
            h_rr_star[i] = max(H_rr_star[i] - b_rr_star, threshold)
        else
            h_ll_star[i] = max(H_ll_star[i] - H_ll_star[i + 1], threshold)
            h_rr_star[i] = max(H_rr_star[i] - H_rr_star[i + 1], threshold)
        end
    end

    # Reconstructed states
    u_ll_star = SVector(h_ll_star..., (h_ll_star .* v1_ll)..., (h_ll_star .* v2_ll)...,
                        b_ll_star)
    u_rr_star = SVector(h_rr_star..., (h_rr_star .* v1_rr)..., (h_rr_star .* v2_rr)...,
                        b_rr_star)

    return SVector(u_ll_star, u_rr_star)
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the bottom
# topography
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::ShallowWaterMultiLayerEquations2D)
    λ = dissipation.max_abs_speed(u_ll, u_rr, orientation_or_normal_direction,
                                  equations)
    diss = -0.5 * λ * (u_rr - u_ll)
    return SVector(@views diss[1:(end - 1)]..., zero(eltype(u_ll)))
end

@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterMultiLayerEquations2D)
    h = waterheight(u, equations)
    h_v1, h_v2 = momentum(u, equations)

    # Calculate averaged velocity of both layers
    H = sum(h)
    v1_m = sum(h_v1) / H
    v2_m = sum(h_v2) / H
    c = sqrt(equations.gravity * H)

    return abs(v1_m) + c, abs(v2_m) + c
end

# Calculate approximation for maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound. This function uses approximate
# eigenvalues as there is no simple way to calculate them analytically.
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr,
                                           orientation::Integer,
                                           equations::ShallowWaterMultiLayerEquations2D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)

    # Get the momentum quantities in the appropriate direction
    if orientation == 1
        h_v_ll, _ = momentum(u_ll, equations)
        h_v_rr, _ = momentum(u_rr, equations)
    else
        _, h_v_ll = momentum(u_ll, equations)
        _, h_v_rr = momentum(u_rr, equations)
    end

    # Get the averaged velocity
    v_m_ll = sum(h_v_ll) / sum(h_ll)
    v_m_rr = sum(h_v_rr) / sum(h_rr)

    # Calculate the wave celerity on the left and right
    c_ll = sqrt(equations.gravity * sum(h_ll))
    c_rr = sqrt(equations.gravity * sum(h_rr))

    return max(abs(v_m_ll), abs(v_m_rr)) + max(c_ll, c_rr)
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr,
                                           normal_direction::AbstractVector,
                                           equations::ShallowWaterMultiLayerEquations2D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    h_v1_ll, h_v2_ll = momentum(u_ll, equations)
    h_v1_rr, h_v2_rr = momentum(u_rr, equations)

    # Get the averaged velocity
    v1_m_ll = sum(h_v1_ll) / sum(h_ll)
    v2_m_ll = sum(h_v2_ll) / sum(h_ll)
    v1_m_rr = sum(h_v1_rr) / sum(h_rr)
    v2_m_rr = sum(h_v2_rr) / sum(h_rr)

    # Compute velocity in the normal direction
    v_dot_n_ll = v1_m_ll * normal_direction[1] + v2_m_ll * normal_direction[2]
    v_dot_n_rr = v1_m_rr * normal_direction[1] + v2_m_rr * normal_direction[2]

    # Calculate the wave celerity on the left and right
    c_ll = sqrt(equations.gravity * sum(h_ll))
    c_rr = sqrt(equations.gravity * sum(h_rr))

    # The normal velocities are already scaled by the norm
    return (max(abs(v_dot_n_ll), abs(v_dot_n_rr)) +
            max(c_ll, c_rr) * norm(normal_direction))
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterMultiLayerEquations2D)
    # Extract waterheight
    h = waterheight(u, equations)
    b = u[end]

    # Initialize total layer height
    H = MVector{nlayers(equations), real(equations)}(undef)
    for i in reverse(eachlayer(equations))
        if i == nlayers(equations)
            setindex!(H, h[i] + b, i)
        else
            setindex!(H, h[i] + H[i + 1], i)
        end
    end

    v1, v2 = velocity(u, equations)
    return SVector{3 * nlayers(equations) + 1, real(equations)}(H..., v1..., v2..., b)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterMultiLayerEquations2D)
    # To extract the total layer height and velocity we reuse the waterheight and momentum functions 
    # from the conservative variables.
    H = waterheight(prim, equations)    # For primitive variables this extracts the total layer height
    v1, v2 = momentum(prim, equations)  # For primitive variables this extracts the velocities
    b = prim[end]

    # Calculate waterheight
    h = MVector{nlayers(equations), real(equations)}(undef)
    for i in eachlayer(equations)
        if i < nlayers(equations)
            h[i] = H[i] - H[i + 1]
        else
            # The lowest layer is measured from the bottom topography
            h[i] = H[i] - b
        end
    end

    # Calculate momentum
    h_v1 = SVector{nlayers(equations), real(equations)}(h[i] * v1[i]
                                                        for i in eachlayer(equations))
    h_v2 = SVector{nlayers(equations), real(equations)}(h[i] * v2[i]
                                                        for i in eachlayer(equations))

    return SVector{3 * nlayers(equations) + 1, real(equations)}(h..., h_v1..., h_v2...,
                                                                b)
end

# Convert conservative variables to entropy variables
# Note, only the first `3 * nlayers(equations)` are the entropy variables
# The last entry still just carries the bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::ShallowWaterMultiLayerEquations2D)
    # Extract conservative variables and compute velocity
    h = waterheight(u, equations)
    b = u[end]
    v1, v2 = velocity(u, equations)
    g = equations.gravity

    # Initialize entropy variable vector
    w = MVector{3 * nlayers(equations) + 1, real(equations)}(undef)

    # Calculate entropy variables in each layer
    for i in eachlayer(equations)
        # Compute w1[i] = ρ[i] * g * (b + ∑h[k] + ∑σ[k] * h[k]) - 0.5 * ρ[i] * (v1[i]^2 + v2[i]^2), 
        # where σ[k] = ρ[k] / ρ[i] denotes the density ratio of different layers
        w1 = equations.rhos[i] * (g * b) - 0.5 * equations.rhos[i] * (v1[i]^2 + v2[i]^2)
        for j in eachlayer(equations)
            if j < i
                w1 += equations.rhos[i] * g *
                      (equations.rhos[j] / equations.rhos[i] * h[j])
            else # i<j<nlayers
                w1 += equations.rhos[i] * g * h[j]
            end
        end

        w2 = equations.rhos[i] * v1[i]
        w3 = equations.rhos[i] * v2[i]

        setlayer!(w, w1, w2, w3, i, equations)
    end
    setindex!(w, b, nlayers(equations) * 3 + 1)
    return SVector(w)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterMultiLayerEquations2D)
    return SVector{nlayers(equations), real(equations)}(u[i]
                                                        for i in 1:nlayers(equations))
end

@inline function momentum(u, equations::ShallowWaterMultiLayerEquations2D)
    h_v1 = SVector{nlayers(equations), real(equations)}(u[i]
                                                        for i in (nlayers(equations) + 1):(2 * nlayers(equations)))
    h_v2 = SVector{nlayers(equations), real(equations)}(u[i]
                                                        for i in (2 * nlayers(equations) + 1):(3 * nlayers(equations)))
    return h_v1, h_v2
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterMultiLayerEquations2D)
    h = waterheight(u, equations)
    h_v1, h_v2 = momentum(u, equations)

    # Compute velocity
    v1 = SVector{nlayers(equations), real(equations)}(h_v1[i] / h[i]
                                                      for i in 1:nlayers(equations))
    v2 = SVector{nlayers(equations), real(equations)}(h_v2[i] / h[i]
                                                      for i in 1:nlayers(equations))
    return v1, v2
end

# Entropy function for the multilayer shallow water equations is the total energy
@inline function Trixi.entropy(u, equations::ShallowWaterMultiLayerEquations2D)
    energy_total(u, equations)
end

# Calculate total energy for a conservative state `u`.
# The total energy is composed of the kinetic energy, the hydrostatic pressure and the potential
# energy from the bottom topography and the layer interfaces.
@inline function Trixi.energy_total(u, equations::ShallowWaterMultiLayerEquations2D)
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)
    b = u[end]
    g = equations.gravity

    e = zero(real(equations))
    for i in eachlayer(equations)
        e += (equations.rhos[i] *
              (0.5 * (h[i] * v1[i]^2 + v2[i]^2 + g * h[i]^2) + g * h[i] * b))
        for j in 1:(i - 1)
            e += g * equations.rhos[j] * h[j] * h[i]
        end
    end

    return e
end

# Calculate kinetic energy for a conservative state `u`
@inline function Trixi.energy_kinetic(u, equations::ShallowWaterMultiLayerEquations2D)
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)

    return (0.5 * sum(equations.rhos[i] * h[i] * (v1[i]^2 + v2[i]^2)
                for i in eachlayer(equations)))
end

# Calculate potential energy for a conservative state `u`
@inline function Trixi.energy_internal(u, equations::ShallowWaterMultiLayerEquations2D)
    return energy_total(u, equations) - energy_kinetic(u, equations)
end

@inline function Trixi.waterheight_pressure(u,
                                            equations::ShallowWaterMultiLayerEquations2D)
    h = waterheight(u, equations)

    return 0.5 * equations.gravity * sum(h .^ 3)
end

# Calculate the error for the "lake-at-rest" test case where H = ∑h + b should
# be a constant value over time. 
# Note, assumes there is a single reference water height `H0` with which to compare.
@inline function Trixi.lake_at_rest_error(u,
                                          equations::ShallowWaterMultiLayerEquations2D)
    h = waterheight(u, equations)
    b = u[end]

    # For well-balancedness testing with possible wet/dry regions the reference
    # water height `H0` accounts for the possibility that the bottom topography
    # can emerge out of the water as well as for the threshold offset to avoid
    # division by a "hard" zero water heights as well.
    H0_wet_dry = max(equations.H0, b + equations.threshold_limiter)

    return abs(H0_wet_dry - (sum(h) + b))
end

# Helper function to set the layer values in the flux computation
@inline function setlayer!(f, f_h, f_hv1, f_hv2, i,
                           equations::ShallowWaterMultiLayerEquations2D)
    setindex!(f, f_h, i)
    setindex!(f, f_hv1, i + nlayers(equations))
    setindex!(f, f_hv2, i + 2 * nlayers(equations))
    return nothing
end
end # @muladd
