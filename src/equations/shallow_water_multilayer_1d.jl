# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    ShallowWaterMultiLayerEquations1D(gravity, H0, rhos)

Multi-Layer Shallow Water equations (MLSWE) in one space dimension. The equations are given by
```math
\left\{
	\begin{aligned}			
		&\partial_t h_m + \partial_x h_mv_m = 0,\\
		&\partial h_mv_m + \partial_x h_mv_m^2 = -gh_m\partial_x \bigg(b + \sum\limits_{k\geq j}h_k + \sum\limits_{k<m}\frac{\rho_k}{\rho_m}h_k \bigg)
	\end{aligned}
\right.
```

where ``m = 1, 2, ..., M`` is the layer index and the unknown variables are the water height ``h`` and
the velocity ``v``.  Furthermore, ``g`` denotes the gravitational constant, ``b(x)`` the bottom 
topography and ``\rho`` the layer density, that must be chosen such that 
``\rho_1 < \rho_2 < ... < \rho_M``, to make sure that layers are ordered from top to bottom, with 
increasing layer density.

The additional quantity ``H_0`` is also available to store a reference value for the total water
height that is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

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

struct ShallowWaterMultiLayerEquations1D{NVARS, NLAYERS, RealT <: Real} <:
       AbstractShallowWaterMultiLayerEquations{1, NVARS, NLAYERS}
    gravity::RealT   # gravitational constant
    H0::RealT        # constant "lake-at-rest" total water height
    rhos::SVector{NLAYERS, RealT} # Vector of layer densities

    function ShallowWaterMultiLayerEquations1D{NVARS, NLAYERS, RealT}(gravity::RealT,
                                                                      H0::RealT,
                                                                      rhos::SVector{NLAYERS,
                                                                                    RealT}) where {
                                                                                                   NVARS,
                                                                                                   NLAYERS,
                                                                                                   RealT <:
                                                                                                   Real
                                                                                                   }
        ## Ensure that layer densities are all positive and in increasing order:
        # Check for increasing order
        issorted(rhos) ||
            throw(ArgumentError("densities must be in increasing order (rhos[1] < rhos[2] < ... < rhos[NLAYERS])"))
        # Check for positive values
        min(rhos...) > 0 || throw(ArgumentError("densities must be positive"))

        new(gravity, H0, rhos)
    end
end

# Allow for flexibility to set the gravitational constant within an elixir depending on the
# application where `gravity_constant=1.0` or `gravity_constant=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases. 
function ShallowWaterMultiLayerEquations1D(; gravity_constant,
                                           H0 = zero(gravity_constant), rhos)

    # Promote all variables to a common type
    _rhos = promote(rhos...)
    RealT = promote_type(eltype(_rhos), eltype(gravity_constant), eltype(H0))

    NLAYERS = length(rhos)
    NVARS = 2 * NLAYERS + 1

    __rhos = SVector(map(RealT, _rhos))
    return ShallowWaterMultiLayerEquations1D{NVARS, NLAYERS, RealT}(gravity_constant,
                                                                    H0, __rhos)
end

@inline function Base.real(::ShallowWaterMultiLayerEquations1D{NVARS, NLAYERS, RealT}) where {
                                                                                              NVARS,
                                                                                              NLAYERS,
                                                                                              RealT <:
                                                                                              Real
                                                                                              }
    RealT
end

Trixi.have_nonconservative_terms(::ShallowWaterMultiLayerEquations1D) = True()
function Trixi.varnames(::typeof(cons2cons),
                        equations::ShallowWaterMultiLayerEquations1D)
    heights = ntuple(n -> "h" * string(n), Val(nlayers(equations)))
    momentas = ntuple(n -> "h" * string(n) * "_v", Val(nlayers(equations)))

    return (heights..., momentas..., "b")
end

# Note, we use the total water heights, H = ∑h + b as primitive variables for easier visualization and setting initial
# conditions
function Trixi.varnames(::typeof(cons2prim),
                        equations::ShallowWaterMultiLayerEquations1D)
    heights = ntuple(n -> "H" * string(n), Val(nlayers(equations)))
    velocities = ntuple(n -> "v" * string(n) * "_1", Val(nlayers(equations)))
    return (heights..., velocities..., "b")
end

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterMultiLayerEquations1D)

A smooth initial condition for a two-layer configuration used for convergence tests in combination with
[`source_terms_convergence_test`](@ref) (and
[`BoundaryConditionDirichlet(initial_condition_convergence_test)`](@ref) in non-periodic domains).
"""
function Trixi.initial_condition_convergence_test(x, t,
                                                  equations::ShallowWaterMultiLayerEquations1D)
    # some constants are chosen such that the function is periodic on the domain [0,sqrt(2)]
    ω = 2.0 * pi * sqrt(2.0)

    H1 = 4.0 + 0.1 * cos(ω * x[1] + t)
    H2 = 2.0 + 0.1 * sin(ω * x[1] + t)
    v1 = 0.9
    v2 = 1.0
    b = 1.0 + 0.1 * cos(2.0 * ω * x[1])

    return prim2cons(SVector(H1, H2, v1, v2, b), equations)
end

"""
    source_terms_convergence_test(u, x, t, equations::ShallowWaterMultiLayerEquations1D)

Source terms used for convergence tests with a two-layer configuration in combination with
[`initial_condition_convergence_test`](@ref)
(and [`BoundaryConditionDirichlet(initial_condition_convergence_test)`](@ref)
in non-periodic domains).
"""
@inline function Trixi.source_terms_convergence_test(u, x, t,
                                                     equations::ShallowWaterMultiLayerEquations1D)
    # Same settings as in `initial_condition_convergence_test`. Some derivative simplify because
    # this manufactured solution velocity is taken to be constant
    ω = 2 * pi * sqrt(2.0)

    du1 = (-0.1 * cos(t + ω * x[1]) - 0.1 * sin(t + ω * x[1]) -
           0.09 * ω * cos(t + ω * x[1]) +
           -0.09 * ω * sin(t + ω * x[1]))
    du2 = 0.1 * cos(t + ω * x[1]) + 0.1 * ω * cos(t + ω * x[1]) +
          0.2 * ω * sin(2.0 * ω * x[1])
    du3 = (5.0 * (-0.1 * ω * cos(t + ω * x[1]) - 0.1 * ω * sin(t + ω * x[1])) *
           (4.0 + 0.2 * cos(t + ω * x[1]) +
            -0.2 * sin(t + ω * x[1])) +
           0.1 * ω * (20.0 + cos(t + ω * x[1]) - sin(t + ω * x[1])) *
           cos(t +
               ω * x[1]) - 0.09 * cos(t + ω * x[1]) - 0.09 * sin(t + ω * x[1]) -
           0.081 * ω * cos(t + ω * x[1]) +
           -0.081 * ω * sin(t + ω * x[1]))
    du4 = ((10.0 + sin(t + ω * x[1]) - cos(2ω * x[1])) *
           (-0.09 * ω * cos(t + ω * x[1]) - 0.09 * ω * sin(t +
                                                           ω * x[1]) -
            0.2 * ω * sin(2 * ω * x[1])) + 0.1 * cos(t + ω * x[1]) +
           0.1 * ω * cos(t + ω * x[1]) +
           5.0 * (0.1 * ω * cos(t + ω * x[1]) + 0.2 * ω * sin(2.0 * ω * x[1])) *
           (2.0 + 0.2 * sin(t + ω * x[1]) +
            -0.2 * cos(2.0 * ω * x[1])) + 0.2 * ω * sin(2.0 * ω * x[1]))

    return SVector(du1, du2, du3, du4, zero(eltype(u)))
end

# Calculate 1D flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::ShallowWaterMultiLayerEquations1D)
    # Extract waterheights and momentas and compute velocities
    hv = momentas(u, equations)
    v = velocity(u, equations)

    # Initialize flux vector
    f = zero(MVector{2 * nlayers(equations) + 1, real(equations)})
    # Calculate fluxes in each layer. 
    # Note that the momentum flux simplifies as the pressure is included in the nonconservative term.
    for i in eachlayer(equations)
        f_h = hv[i]
        f_hv = hv[i] * v[i]

        setlayer!(f, f_h, f_hv, i, equations)
    end
    return SVector(f)
end

"""
    flux_nonconservative_ersing_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterTwoLayerEquations1D)

Non-symmetric path-conservative two-point flux discretizing the nonconservative (source) term
that contains the gradients of the bottom topography and waterheights from the coupling between layers
and the nonconservative pressure formulation [`ShallowWaterMultiLayerEquations1D`](@ref).

When the bottom topography is nonzero this scheme will be well-balanced when used with the 
nonconservative [`flux_ersing_etal`](@ref).

In the two-layer setting this combination is equivalent to the fluxes in:
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on 
  curvilinear meshes
  [DOI: 10.48550/arXiv.2306.12699](https://doi.org/10.48550/arXiv.2306.12699)
"""
@inline function Trixi.flux_nonconservative_ersing_etal(u_ll, u_rr,
                                                        orientation::Integer,
                                                        equations::ShallowWaterMultiLayerEquations1D)
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
    f = zero(MVector{2 * nlayers(equations) + 1, real(equations)})

    # Compute the nonconservative flux in each layer (0, ..., 0, f_hv[1], ..., f_hv[NLAYERS], 0)
    # with f_hv[i] = gh[i] * (b + ∑h[k] + ∑σ[k]h[k])_x.  
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
        setindex!(f, f_hv, i + nlayers(equations))
    end

    return SVector(f)
end

"""
    flux_ersing_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterMultiLayerEquations1D)

Total energy conservative (mathematical entropy for MLSWE) split form,
without the hydrostatic pressure.
When the bottom topography is nonzero this scheme will be well-balanced when used with the 
nonconservative [`flux_nonconservative_ersing_etal`](@ref).

In the two-layer setting this combination is equivalent to the fluxes in:
- Patrick Ersing, Andrew R. Winters (2023)
  An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on 
  curvilinear meshes
  [DOI: 10.48550/arXiv.2306.12699](https://doi.org/10.48550/arXiv.2306.12699)
"""
@inline function flux_ersing_etal(u_ll, u_rr,
                                  orientation::Integer,
                                  equations::ShallowWaterMultiLayerEquations1D)
    # Unpack left and right state
    hv_ll = momentas(u_ll, equations)
    hv_rr = momentas(u_rr, equations)

    # Get the velocities on either side
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Initialize flux vector
    f = zero(MVector{2 * nlayers(equations) + 1, real(equations)})

    # Calculate fluxes in each layer. 
    # Note that the momentum flux simplifies as the pressure is included in the nonconservative term.
    for i in eachlayer(equations)
        # Compute averages
        v_avg = 0.5 * (v_ll[i] + v_rr[i])
        hv_avg = 0.5 * (hv_ll[i] + hv_rr[i])

        f_h = hv_avg
        f_hv = f_h * v_avg

        setlayer!(f, f_h, f_hv, i, equations)
    end

    return SVector(f)
end

@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterMultiLayerEquations1D)
    h = waterheight(u, equations)
    hv = momentas(u, equations)

    # Calculate averaged velocity of both layers
    H = sum(h)
    v_m = sum(hv) / H
    c = sqrt(equations.gravity * H)

    return (abs(v_m) + c)
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterMultiLayerEquations1D)
    # Extract waterheights and momentas
    h = waterheight(u, equations)
    b = u[end]

    # Initialize total waterheight
    H = MVector{nlayers(equations), real(equations)}(undef)
    for i in reverse(eachlayer(equations))
        if i == nlayers(equations)
            setindex!(H, h[i] + b, i)
        else
            setindex!(H, h[i] + H[i + 1], i)
        end
    end

    v = velocity(u, equations)
    return SVector{2 * nlayers(equations) + 1, real(equations)}(H..., v..., b)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterMultiLayerEquations1D)
    H = prim[1:nlayers(equations)]
    v = prim[(nlayers(equations) + 1):(2 * nlayers(equations))]
    b = prim[end]

    # Calculate waterheights
    h = MVector{nlayers(equations), real(equations)}(undef)
    for i in eachlayer(equations)
        if i < nlayers(equations)
            setindex!(h, H[i] - H[i + 1], i)
        else
            # The lowest layer is measured from the bottom topography
            setindex!(h, H[i] - b, i)
        end
    end

    # Calculate momentas
    h_v = SVector{nlayers(equations), real(equations)}(h[i] * v[i]
                                                       for i in eachlayer(equations))

    return SVector{2 * nlayers(equations) + 1, real(equations)}(h..., h_v..., b)
end

# Convert conservative variables to entropy variables
# Note, only the first four are the entropy variables, the fifth entry still just carries the
# bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::ShallowWaterMultiLayerEquations1D)
    # Extract conservative variables and compute velocity
    h = waterheight(u, equations)
    b = u[end]
    v = velocity(u, equations)
    g = equations.gravity

    # Initialize entropy variable vector
    w = MVector{2 * nlayers(equations) + 1, real(equations)}(undef)

    # Calculate entropy variables in each layer
    for i in eachlayer(equations)
        # Compute w1[i] = ρ[i]g * (b + ∑h[k] + ∑σ[k]h[k])
        w1 = equations.rhos[i] * g * b
        for j in eachlayer(equations)
            if j < i
                w1 += equations.rhos[i] * g *
                      (equations.rhos[j] / equations.rhos[i] * h[j])
            else # i<j<nlayers
                w1 += equations.rhos[i] * g * h[j]
            end
        end

        w2 = equations.rhos[i] * v[i]

        setlayer!(w, w1, w2, i, equations)
    end
    setindex!(w, b, nlayers(equations) * 2 + 1)
    return SVector(w)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterMultiLayerEquations1D)
    return SVector{nlayers(equations), real(equations)}(u[i]
                                                        for i in 1:nlayers(equations))
end

@inline function momentas(u, equations::ShallowWaterMultiLayerEquations1D)
    return SVector{nlayers(equations), real(equations)}(u[i]
                                                        for i in (nlayers(equations) + 1):(2 * nlayers(equations)))
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterMultiLayerEquations1D)
    h = waterheight(u, equations)
    hv = momentas(u, equations)
    return SVector{nlayers(equations), real(equations)}(hv[i] / h[i]
                                                        for i in 1:nlayers(equations))
end

# Helper function to set the layer values in the flux computation
@inline function setlayer!(f, f_h, f_hv, i,
                           equations::ShallowWaterMultiLayerEquations1D)
    setindex!(f, f_h, i)
    setindex!(f, f_hv, i + nlayers(equations))
end

# Calculate the error for the "lake-at-rest" test case where H = ∑h+b should
# be a constant value over time
@inline function Trixi.lake_at_rest_error(u,
                                          equations::ShallowWaterMultiLayerEquations1D)
    h = waterheight(u, equations)
    b = u[end]

    return abs(equations.H0 - (sum(h) + b))
end
end # @muladd
