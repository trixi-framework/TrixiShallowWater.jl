# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    ShallowWaterMomentEquations1D(; gravity, H0 = zero(gravity), n_moments,
                                    nu = 0.1, lambda = 0.1, rho = 1000.0,
                                    nman = 0.0165)
Shallow water moment equations in one spatial dimensions. The equations are given by
```math
\begin{cases}
\partial_t h + \partial_x hv = 0, \\
\partial_t hv + \partial_x (hv^2 + h\sum\limits_{i=1}^N \frac{\alpha_i^2}{2i+1}) = -gh\partial(h+b)_x,\\
\partial_t h\alpha_i + \partial_x \left(2hv\alpha_i + h\sum\limits_{j,k=1}^N A_{ijk}\alpha_j \alpha_k \right) = 
                       v\partial_x h\alpha_i - h\sum\limits_{j,k=1}^N B_{ijk} \alpha_k \partial_x h\alpha_j,
\end{cases}
```

The unknown quantities are the water and sediment height ``h``, the velocity ``v`` and the moments
``\alpha_i`` for ``i = 1, ..., n_moments``. The terms ``A_{ijk}`` and ``B_{ijk}`` are moment tensors that 
are precomputed using shifted Legendre polynomials and ``g`` is the gravitational acceleration.

The conservative variable water height ``h`` is measured from the bottom topography ``b``, therefore
one also defines the total water height as ``H = h + b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water
height that is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

Additionally the system can be extended with two different friction models.
A Newtonian slip friction model is available that uses slip length ``\lambda`` and kinematic viscosity ``\nu`` as parameters.
The source term for this friction model is given by [`source_term_newtonian_slip_friction`](@ref).
A Manning friction model is also available that uses the Manning roughness coefficient `nman` and
fluid density `rho` as a parameter.
The source term for this friction model is given by [`source_term_manning_friction`](@ref).

For details see the paper:
- Julio Careaga, Patrick Ersing, Julian Koellermeier, Andrew R. Winters (2026)
  Entropy analysis and entropy stable DG methods for the shallow water moment equations
  [DOI: 10.48550/arXiv.2602.06513](https://doi.org/10.48550/arXiv.2602.06513)
"""
struct ShallowWaterMomentEquations1D{NVARS, NMOMENTS, RealT <: Real} <:
       AbstractMomentEquations{1, NVARS, NMOMENTS}
    gravity::RealT   # gravitational acceleration
    H0::RealT        # constant "lake-at-rest" total water height
    n_moments::Integer  # number of moments
    # Moment matrices
    A::Array{RealT, 3}
    B::Array{RealT, 3}
    C::Array{RealT, 2}
    # Friction related quantities
    nu::RealT       # kinematic viscosity
    lambda::RealT  # slip length
    rho::RealT    # fluid density (relevant for Manning friction)
    nman::RealT   # Manning roughness coefficient

    function ShallowWaterMomentEquations1D{NVARS, NMOMENTS, RealT}(gravity::RealT,
                                                                   H0::RealT,
                                                                   n_moments::Integer,
                                                                   A,
                                                                   B::Array{RealT, 3},
                                                                   C::Array{RealT, 2},
                                                                   nu::RealT,
                                                                   lambda::RealT,
                                                                   rho::RealT,
                                                                   nman::RealT) where {
                                                                                       NVARS,
                                                                                       NMOMENTS,
                                                                                       RealT <:
                                                                                       Real
                                                                                       }
        new(gravity, H0, n_moments, A, B, C, nu, lambda, rho, nman)
    end
end

# Allow for flexibility to set the gravitational acceleration and number of moments within an elixir
# depending on the application. Here `gravity=1.0` or `gravity=9.81` are common values for the
# gravitational acceleration. The reference total water height H0 defaults to 0.0 but is used for 
# the "lake-at-rest" well-balancedness test cases.
function ShallowWaterMomentEquations1D(; gravity, H0 = zero(gravity), n_moments,
                                       nu = 0.1, lambda = 0.1, rho = 1000.0,
                                       nman = 0.0165)
    RealT = promote_type(typeof(gravity), typeof(H0))

    # Extract number of moments and variables
    NMOMENTS = n_moments
    NVARS = NMOMENTS + 3    # (h, hv, a_1, ..., a_NMOMENTS, b)

    # Compute moment matrices
    B = compute_B_tensor(n_moments)
    A = compute_A_tensor(B)
    C = compute_C_matrix(n_moments)

    return ShallowWaterMomentEquations1D{NVARS, NMOMENTS, RealT}(gravity,
                                                                 H0,
                                                                 n_moments,
                                                                 A,
                                                                 B,
                                                                 C,
                                                                 nu,
                                                                 lambda,
                                                                 rho,
                                                                 nman)
end

@inline function Base.real(::ShallowWaterMomentEquations1D{NVARS, NMOMENTS, RealT}) where {
                                                                                           NVARS,
                                                                                           NMOMENTS,
                                                                                           RealT <:
                                                                                           Real
                                                                                           }
    RealT
end

Trixi.have_nonconservative_terms(::ShallowWaterMomentEquations1D) = True()

function Trixi.varnames(::typeof(cons2cons), equations::ShallowWaterMomentEquations1D)
    conservative_moments = ntuple(n -> "ha" * string(n), Val(nmoments(equations)))

    return ("h", "hv", conservative_moments..., "b")
end

# We use the total layer heights, H = h + b as a primitive variables for easier visualization and setting initial
# conditions
function Trixi.varnames(::typeof(cons2prim), equations::ShallowWaterMomentEquations1D)
    primitive_moments = ntuple(n -> "a" * string(n), Val(nmoments(equations)))
    return ("H", "v", primitive_moments..., "b")
end

"""
    source_term_newtonian_slip_friction(u, x, t, equations::Union{ShallowWaterMomentEquations1D,
                                                    ShallowWaterLinearizedMomentEquations1D})
  
  Bottom friction source term using a Newtonian slip friction model. 
  The parameters slip length `lambda` and kinematic viscosity `nu` are set in the equation struct.

  For details see the paper:
  - Julio Careaga, Patrick Ersing, Julian Koellermeier, Andrew R. Winters (2026)
    Entropy analysis and entropy stable DG methods for the shallow water moment equations
    [DOI: 10.48550/arXiv.2602.06513](https://doi.org/10.48550/arXiv.2602.06513)
"""
@inline function source_term_newtonian_slip_friction(u, x, t,
                                                     equations::Union{ShallowWaterMomentEquations1D,
                                                                      ShallowWaterLinearizedMomentEquations1D})
    # Get waterheight, velocity and moments
    h = waterheight(u, equations)
    v = velocity(u, equations)
    a = moments(u, equations) / h

    sum_a = sum(a)

    friction_v = -equations.nu / equations.lambda * (v + sum_a)

    friction_mom = MVector{nmoments(equations), real(equations)}(undef)

    for i in eachmoment(equations)
        friction_mom[i] = -(2 * i + 1) * equations.nu / equations.lambda * (v + sum_a)
        for j in eachmoment(equations)
            friction_mom[i] += -equations.nu / h * equations.C[i, j] * a[j]
        end
    end

    return SVector(0, friction_v, friction_mom..., 0)
end

"""
    source_term_manning_friction(u, x, t, equations::Union{ShallowWaterMomentEquations1D,
                                                    ShallowWaterLinearizedMomentEquations1D})

  Bottom friction source term using a Manning friction model.
  The Manning roughness coefficient `nman` is set in the equation struct.

  For details see the paper:
  - Julio Careaga, Patrick Ersing, Julian Koellermeier, Andrew R. Winters (2026)
    Entropy analysis and entropy stable DG methods for the shallow water moment equations
    [DOI: 10.48550/arXiv.2602.06513](https://doi.org/10.48550/arXiv.2602.06513)
"""
@inline function source_term_manning_friction(u, x, t,
                                              equations::Union{ShallowWaterMomentEquations1D,
                                                               ShallowWaterLinearizedMomentEquations1D})
    # Get waterheight, velocity and moments
    h = waterheight(u, equations)
    v = velocity(u, equations)
    a = moments(u, equations) / h

    sum_a = sum(a)

    friction_v = -equations.rho * equations.gravity * ((equations.nman^2) / h^(1 / 3)) *
                 (v + sum_a) * abs(v + sum_a)

    friction_mom = MVector{nmoments(equations), real(equations)}(undef)

    for i in eachmoment(equations)
        friction_mom[i] = -(2 * i + 1) * equations.rho * equations.gravity *
                          ((equations.nman^2) / h^(1 / 3)) * (v + sum_a) *
                          abs(v + sum_a)
        for j in eachmoment(equations)
            friction_mom[i] += -equations.nu / h * equations.C[i, j] * a[j]
        end
    end

    return SVector(0, friction_v, friction_mom..., 0)
end

"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, x, t, surface_flux_function,
                                 equations::ShallowWaterMomentEquations1D)

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
                                                    orientation_or_normal,
                                                    direction,
                                                    x,
                                                    t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterMomentEquations1D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # Create the "external" boundary solution state
    h = u_inner[1]
    hv = u_inner[2]
    ha = moments(u_inner, equations)
    b = u_inner[end]

    u_boundary = SVector(h, -hv, ha..., b)

    # Calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal,
                                     equations)
        noncons_flux = nonconservative_flux_function(u_inner,
                                                     u_boundary,
                                                     orientation_or_normal,
                                                     equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal,
                                     equations)
        noncons_flux = nonconservative_flux_function(u_inner,
                                                     u_boundary,
                                                     orientation_or_normal,
                                                     equations)
    end
    return flux, noncons_flux
end

# Calculate 1D advective portion of the flux for a single point
@inline function Trixi.flux(u,
                            orientation::Integer,
                            equations::ShallowWaterMomentEquations1D)
    # Extract conservative variables    
    h = waterheight(u, equations)
    hv = u[2]
    ha = moments(u, equations)

    # Compute primitive variables
    v = velocity(u, equations)
    a = ha / h

    f1 = hv
    f2 = hv * v
    for i in eachmoment(equations)
        f2 += ha[i] * a[i] / (2 * i + 1)
    end

    # Compute the moment fluxes
    f_moments = MArray{nmoments(equations), real(equations)}(2 * ha * v)

    # Add the moment tensor contributions
    for i in eachmoment(equations),
        j in eachmoment(equations),
        k in eachmoment(equations)

        f_moments[i] += equations.A[i, j, k] * a[j] * a[k]
    end

    return SVector{nmoments(equations) + 3, real(equations)}(f1, f2, f_moments..., 0)
end

"""
    flux_careaga_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterMomentEquations1D)

Total energy conservative split form, without the hydrostatic pressure.
When the bottom topography is nonzero this scheme will be well-balanced when used with the 
nonconservative flux [`flux_nonconservative_careaga_etal`](@ref).

To obtain an entropy stable and well-balanced formulation the `surface_flux` can be set as
`FluxPlusDissipation(flux_careaga_etal, DissipationLaxFriedrichsEntropyVariables()), flux_nonconservative_careaga_etal`.

For details see the paper:
  - Julio Careaga, Patrick Ersing, Julian Koellermeier, Andrew R. Winters (2026)
    Entropy analysis and entropy stable DG methods for the shallow water moment equations
    [DOI: 10.48550/arXiv.2602.06513](https://doi.org/10.48550/arXiv.2602.06513)
"""
@inline function flux_careaga_etal(u_ll,
                                   u_rr,
                                   orientation::Integer,
                                   equations::ShallowWaterMomentEquations1D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    hv_ll = u_ll[2]
    hv_rr = u_rr[2]
    ha_ll = moments(u_ll, equations)
    ha_rr = moments(u_rr, equations)

    # Get the velocities and primitive moments on either side
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)
    a_ll = ha_ll / h_ll
    a_rr = ha_rr / h_rr

    v_avg = 0.5 * (v_ll + v_rr)
    hv_avg = 0.5 * (hv_ll + hv_rr)
    a_avg = 0.5 * (a_ll + a_rr)
    ha_avg = 0.5 * (ha_ll + ha_rr)

    # Compute the flux components
    f1 = hv_avg
    f2 = hv_avg * v_avg
    for i in eachmoment(equations)
        f2 += ha_avg[i] * a_avg[i] / (2 * i + 1)
    end

    # Compute the moment fluxes
    f_moments = MVector{nmoments(equations), real(equations)}(hv_avg * a_avg +
                                                              ha_avg * v_avg)
    # Add the moment tensor contributions
    for i in eachmoment(equations),
        j in eachmoment(equations),
        k in eachmoment(equations)

        f_moments[i] += equations.A[i, j, k] * ha_avg[j] * a_avg[k]
    end

    return SVector{nmoments(equations) + 3, real(equations)}(f1, f2, f_moments..., 0)
end

"""
    flux_nonconservative_careaga_etal(u_ll, u_rr, orientation::Integer,
                            equations::ShallowWaterMomentEquations1D)

Non-symmetric path-conservative two-point flux discretizing the nonconservative term
that contains the gradients of the bottom topography, the nonconservative pressure formulation and 
the moment tensor contributions.
When the bottom topography is nonzero this scheme will be well-balanced when used with [`flux_careaga_etal`](@ref).

For details see the paper:
  - Julio Careaga, Patrick Ersing, Julian Koellermeier, Andrew R. Winters (2026)
    Entropy analysis and entropy stable DG methods for the shallow water moment equations
    [DOI: 10.48550/arXiv.2602.06513](https://doi.org/10.48550/arXiv.2602.06513)
"""
@inline function flux_nonconservative_careaga_etal(u_ll,
                                                   u_rr,
                                                   orientation::Integer,
                                                   equations::ShallowWaterMomentEquations1D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    ha_ll = moments(u_ll, equations)
    ha_rr = moments(u_rr, equations)
    a_ll = ha_ll / h_ll
    b_ll = u_ll[end]
    b_rr = u_rr[end]

    # Get the velocities and primitive moments on either side
    v_ll = velocity(u_ll, equations)

    # Compute the jumps
    h_jump = h_rr - h_ll
    b_jump = b_rr - b_ll
    ha_jump = ha_rr - ha_ll

    # Get the gravitational acceleration
    g = equations.gravity

    # Compute nonconservative flux components
    f2 = g * h_ll * (h_jump + b_jump)

    # Compute the moment fluxes.
    f_moments = MVector{nmoments(equations), real(equations)}(-v_ll * ha_jump)

    # Add the moment tensor contributions
    for i in eachmoment(equations),
        j in eachmoment(equations),
        k in eachmoment(equations)

        f_moments[i] += equations.B[i, j, k] * a_ll[k] * ha_jump[j]
    end

    return SVector{nmoments(equations) + 3, real(equations)}(0, f2, f_moments..., 0)
end

# Specialized `DissipationLocalLaxFriedrichs` to avoid spurious dissipation in the bottom
# topography. For nonzero bottom topography [`Trixi.DissipationLaxFriedrichsEntropyVariables`](@extref)
# should be used instead as this version is not well-balanced.
@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll,
                                                              u_rr,
                                                              orientation_or_normal_direction,
                                                              equations::ShallowWaterMomentEquations1D)
    λ = dissipation.max_abs_speed(u_ll,
                                  u_rr,
                                  orientation_or_normal_direction,
                                  equations)
    diss = -0.5 * λ * (u_rr - u_ll)
    return SVector(@views diss[1:(end - 1)]..., zero(eltype(u_ll)))
end

# Provably entropy stable and well-balanced local Lax-Friedrichs dissipation for the SWME that avoids
# spurious dissipation in the bottom topography.
# 
# For details see the paper:
#   - Julio Careaga, Patrick Ersing, Julian Koellermeier, Andrew R. Winters (2026)
#     Entropy analysis and entropy stable DG methods for the shallow water moment equations
#     [DOI: 10.48550/arXiv.2602.06513](https://doi.org/10.48550/arXiv.2602.06513) 
@inline function (dissipation::DissipationLaxFriedrichsEntropyVariables)(u_ll,
                                                                         u_rr,
                                                                         orientation_or_normal_direction,
                                                                         equations::Union{ShallowWaterLinearizedMomentEquations1D,
                                                                                          ShallowWaterMomentEquations1D})
    λ = dissipation.max_abs_speed(u_ll,
                                  u_rr,
                                  orientation_or_normal_direction,
                                  equations)

    # Convert to entropy variables
    w_ll = Trixi.cons2entropy(u_ll, equations)
    w_rr = Trixi.cons2entropy(u_rr, equations)
    # Exclude the bottom topography from the jump
    w_jump = SVector{nmoments(equations) + 2}(@views (w_rr - w_ll)[1:(end - 1)])

    # Compute the matrix H = du/dw at the average state
    h_avg = 0.5 * (u_ll[1] + u_rr[1])
    v_avg = 0.5 * (Trixi.velocity(u_ll, equations) + Trixi.velocity(u_rr, equations))
    a_avg = SVector(0.5 * (moments(u_ll, equations) / u_ll[1] +
                     moments(u_rr, equations) / u_rr[1]))
    g = equations.gravity

    # Construct the H matrix from H = yy' + Iz
    y = SVector{nmoments(equations) + 2, real(equations)}(1, v_avg, a_avg...)
    z = SVector{nmoments(equations) + 2, real(equations)}(ntuple(k -> k == 1 ? 0 :
                                                                      k == 2 ?
                                                                      g * h_avg :
                                                                      (2k - 3) * g *
                                                                      h_avg,
                                                                 nmoments(equations) +
                                                                 2))

    H = SMatrix{nmoments(equations) + 2, nmoments(equations) + 2, real(equations)}(1 /
                                                                                   g *
                                                                                   (y *
                                                                                    y' +
                                                                                    diagm(z)))

    diss = SVector{nmoments(equations) + 2, real(equations)}(-0.5 .* λ .* (H * w_jump))

    return SVector{nmoments(equations) + 3, real(equations)}(diss...,
                                                             zero(eltype(u_ll)))
end

# The eigenvalues are approximate with those of the SWLME taken from:
# - Julian Koellermeier, Ernesto Pimentel-García (2022)
#   Steady states and well-balanced schemes for shallow water moment equations with topography
#   [DOI: 10.1016/j.amc.2022.127166](https://doi.org/10.1016/j.amc.2022.127166)
@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterMomentEquations1D)
    h = waterheight(u, equations)
    ha = moments(u, equations)
    v_m = velocity(u, equations)

    # calculate the wave celerity
    c = equations.gravity * h
    for i in eachmoment(equations)
        c += (3 * ha[i] / h) / (2 * i + 1)
    end
    c = sqrt(c)

    return (abs(v_m) + c)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function Trixi.max_abs_speed(u_ll,
                                     u_rr,
                                     orientation::Integer,
                                     equations::ShallowWaterMomentEquations1D)
    # Unpack left and right state
    h_ll = waterheight(u_ll, equations)
    h_rr = waterheight(u_rr, equations)
    ha_ll = moments(u_ll, equations)
    ha_rr = moments(u_rr, equations)
    v_m_rr = velocity(u_rr, equations)
    v_m_ll = velocity(u_ll, equations)

    # Calculate the wave celerity on the left and right
    c_ll = equations.gravity * h_ll
    c_rr = equations.gravity * h_rr
    for i in eachmoment(equations)
        c_ll += (3 * ha_ll[i] / h_ll) / (2 * i + 1)
        c_rr += (3 * ha_rr[i] / h_rr) / (2 * i + 1)
    end
    c_ll = sqrt(c_ll)
    c_rr = sqrt(c_rr)

    return (max(abs(v_m_ll) + c_ll, abs(v_m_rr) + c_rr))
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterMomentEquations1D)
    # Extract conservative variables
    h = waterheight(u, equations)
    ha = moments(u, equations)
    b = u[end]

    # Compute the total water height, velocity and primitive moments
    H = h + b
    a = ha / h
    v = velocity(u, equations)

    return SVector{nmoments(equations) + 3, real(equations)}(H, v, a..., b)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterMomentEquations1D)
    # To extract the total layer height and velocity we reuse the water height and momentum functions 
    # from the conservative variables.
    H = waterheight(prim, equations)
    v = prim[2]
    a = moments(prim, equations)     # For primitive variables this extracts the primitive moments
    b = prim[end]

    # Compute the conservative variables
    h = H - b
    ha = h * a
    hv = h * v

    return SVector{nmoments(equations) + 3, real(equations)}(h, hv, ha..., b)
end

# Convert conservative variables to entropy variables
@inline function Trixi.cons2entropy(u, equations::ShallowWaterMomentEquations1D)
    # Extract conservative variables and compute velocity
    h = waterheight(u, equations)
    v = velocity(u, equations)
    a = moments(u, equations) / h
    b = u[end]
    g = equations.gravity

    # Calculate entropy variables
    w1 = g * (h + b) - 0.5 * v^2
    # add moment contributions
    for i in eachmoment(equations)
        w1 -= 0.5 * a[i]^2 / (2 * i + 1)
    end

    w2 = v

    w_moments = MVector{nmoments(equations), real(equations)}(undef)
    for i in eachmoment(equations)
        w_moments[i] = a[i] / (2 * i + 1)
    end

    return SVector(w1, w2, w_moments..., b)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterMomentEquations1D)
    return u[1]
end

@inline function moments(u, equations::ShallowWaterMomentEquations1D)
    return SVector{nmoments(equations), real(equations)}(u[i]
                                                         for i in (3:(nmoments(equations) + 2)))
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterMomentEquations1D)
    h = waterheight(u, equations)
    hv = u[2]

    return hv / h
end

# The entropy function for the shallow water moment equations is the total energy
@inline function Trixi.entropy(u, equations::ShallowWaterMomentEquations1D)
    h = waterheight(u, equations)
    v = velocity(u, equations)
    a = moments(u, equations) ./ h
    b = u[end]
    g = equations.gravity

    # Calculate the total energy for the mean flow
    e = 0.5 * (h * v^2 + g * h^2) + g * h * b

    # add the contribution from moments
    for i in eachmoment(equations)
        e += 0.5 * h * a[i]^2 / (2 * i + 1)
    end

    return e
end

# Calculate the error for the "lake-at-rest" test case where H = h + b should
# be a constant value over time. 
# Note, assumes there is a single reference water height `H0` with which to compare.
@inline function Trixi.lake_at_rest_error(u, equations::ShallowWaterMomentEquations1D)
    h = waterheight(u, equations)
    b = u[end]

    return abs(equations.H0 - (h + b))
end
end # @muladd
