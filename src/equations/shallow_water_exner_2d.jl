# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    ShallowWaterExnerEquations2D(;gravity, H0 = 0.0,
                                 friction = ManningFriction(n = 0.0),
                                 sediment_model,
                                 porosity,
                                 rho_f, rho_s)
Formulation of the Shallow water-Exner equations in two space dimensions that possesses a mathematical
entropy inequality. The equations are given by
```math
\begin{cases}
\partial_t h + \partial_x hv_1 + \partial_y hv_2 = 0, \\
\partial_t hv_1 + \partial_x (hv_1^2) + \partial_y (hv_1v_2) + gh\partial_x (h + h_b) + g\frac{1}{r}h_s\partial_x (rh + h_b) + \frac{\tau_1}{\rho_f} = 0,\\
\partial_t hv_2 + \partial_x (hv_1v_2) + \partial_y (hv_2^2) + gh\partial_y (h + h_b) + g\frac{1}{r}h_s\partial_y (rh + h_b) + \frac{\tau_2}{\rho_f} = 0,\\
\partial_t h_b + \partial_x q_{s,1} + \partial_y q_{s,2} = 0,
\end{cases}
```
The unknown quantities are the water and sediment height ``h``, ``h_b`` and the velocities ``v_1``, ``v_2``.
The sediment discharges ``q_{s,1}(h, hv)`` and ``q_{s,2}(h, hv)`` are determined by the `sediment_model`.
The sediment model also determines the active sediment height, e.g., ``h_s = q_s1 / v_1``.
Furthermore ``\tau`` denotes the shear stress at the water-sediment interface and is determined by
the `friction` model.
The gravitational acceleration is denoted by ``g``, and ``\rho_f`` and ``\rho_s`` are the fluid and sediment
densities, respectively. The density ratio is given by ``r = \rho_f / \rho_s``, where ``r`` lies between
``0 < r < 1`` as the fluid density ``\rho_f`` should be smaller than the sediment density ``\rho_s``.

The conservative variable water height ``h`` is measured from the sediment height ``h_b``, therefore
one also defines the total water height as ``H = h + h_b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water
height that is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

The entropy conservative formulation in 1D has been derived in the paper:
- E.D. Fernández-Nieto, T.M. de Luna, G. Narbona-Reina and J. de Dieu Zabsonré (2017)\
  Formal deduction of the Saint-Venant–Exner model including arbitrarily sloping sediment beds and
  associated energy\
  [DOI: 10.1051/m2an/2016018](https://doi.org/10.1051/m2an/2016018)
"""
struct ShallowWaterExnerEquations2D{RealT <: Real,
                                    FrictionT <: Friction{RealT},
                                    SedimentT <: SedimentModel{RealT}} <:
       Trixi.AbstractShallowWaterEquations{2, 4}
    # This structure ensures that friction and sediment models are concrete types
    # to prevent allocations
    gravity::RealT # gravitational acceleration
    H0::RealT      # constant "lake-at-rest" total water height
    friction::FrictionT
    sediment_model::SedimentT
    porosity_inv::RealT  # 1/(1-porosity)
    rho_f::RealT    # density of fluid
    rho_s::RealT    # density of sediment
    r::RealT       # density ratio
end

# Allow for flexibility to set the gravitational acceleration within an elixir depending on the
# application where `gravity=1.0` or `gravity=9.81` are common values.
# The reference total water height H0 defaults to 0.0 but is used for the "lake-at-rest"
# well-balancedness test cases.
function ShallowWaterExnerEquations2D(; gravity, H0 = zero(gravity),
                                      friction = ManningFriction(n = 0.0),
                                      sediment_model,
                                      porosity, rho_f, rho_s)
    # Precompute common expressions for the porosity and density ratio
    porosity_inv = inv(1 - porosity)
    r = rho_f / rho_s
    return ShallowWaterExnerEquations2D(gravity, H0, friction, sediment_model,
                                        porosity_inv, rho_f, rho_s, r)
end

Trixi.have_nonconservative_terms(::ShallowWaterExnerEquations2D) = True()
Trixi.varnames(::typeof(cons2cons), ::ShallowWaterExnerEquations2D) = ("h", "hv1",
                                                                       "hv2", "h_b")
# Note, we use the total water height, H = h + h_b, as the first primitive variable for easier
# visualization and setting initial conditions
Trixi.varnames(::typeof(cons2prim), ::ShallowWaterExnerEquations2D) = ("H", "v1", "v2",
                                                                       "h_b")

"""
    boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
                                 surface_flux_function, equations::ShallowWaterExnerEquations2D)

Should be used together with [`Trixi.TreeMesh`](@extref).
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation,
                                                    direction, x, t,
                                                    surface_flux_functions,
                                                    equations::ShallowWaterExnerEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # get the appropriate normal vector from the orientation
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
    source_term_bottom_friction(u, x, t, equations::ShallowWaterExnerEquations2D)

Source term that accounts for the bottom friction in the [ShallowWaterExnerEquations2D](@ref).
The actual friction law is determined through the friction model in `equations.friction`.
"""
@inline function source_term_bottom_friction(u, x, t,
                                             equations::ShallowWaterExnerEquations2D)
    tau1, tau2 = shear_stress(u, equations)
    return SVector(zero(eltype(u)),
                   -tau1,
                   -tau2,
                   zero(eltype(u)))
end

# Calculate 1D flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::ShallowWaterExnerEquations2D)
    _, hv1, hv2, _ = u
    v1, v2 = velocity(u, equations)
    q_s1, q_s2 = sediment_discharge(u, equations)

    if orientation == 1
        f1 = hv1
        f2 = hv1 * v1
        f3 = hv2 * v1
        f4 = q_s1
    else # orientation == 2
        f1 = hv2
        f2 = hv1 * v2
        f3 = hv2 * v2
        f4 = q_s2
    end

    return SVector(f1, f2, f3, f4)
end

"""
    flux_nonconservative_ersing_etal(u_ll, u_rr, orientation, equations::ShallowWaterExnerEquations2D)

Non-symmetric path-conservative two-point flux discretizing the nonconservative terms of the
[`ShallowWaterExnerEquations2D`](@ref) which consists of the hydrostatic pressure of the fluid
layer and an additional pressure contribution from the sediment layer to obtain an entropy
inequality.

This non-conservative flux should be used together with [`flux_ersing_etal`](@ref) to create a
scheme that is entropy conservative and well-balanced.
"""
@inline function flux_nonconservative_ersing_etal(u_ll, u_rr,
                                                  orientation::Integer,
                                                  equations::ShallowWaterExnerEquations2D)
    # Pull the necessary left and right state information
    h_ll, _, _, h_b_ll = u_ll
    h_rr, _, _, h_b_rr = u_rr

    # Calculate jumps
    h_jump = h_rr - h_ll
    h_b_jump = h_b_rr - h_b_ll

    # Compute the effective sediment height
    h_s_ll = effective_sediment_height(u_ll, equations)

    z = zero(eltype(u_ll))

    f = equations.gravity * h_ll * (h_jump + h_b_jump)
    # Additional nonconservative term to obtain entropy conservative formulation
    f += (equations.gravity / equations.r * h_s_ll *
          (equations.r * h_jump + h_b_jump))

    if orientation == 1
        return SVector(z, f, z, z)
    else # orientation == 2
        return SVector(z, z, f, z)
    end
end

"""
    flux_ersing_etal(u_ll, u_rr, orientation::Integer,
                                     equations::ShallowWaterExnerEquations2D)

Entropy conservative split form, without the hydrostatic pressure. This flux should be used
together with the nonconservative [`flux_nonconservative_ersing_etal`](@ref) to create a scheme
that is entropy conservative and well-balanced.

To obtain an entropy stable formulation the `surface_flux` can be set as
`FluxPlusDissipation(flux_ersing_etal, DissipationLocalLaxFriedrichs()), flux_nonconservative_ersing_etal`.
"""
@inline function flux_ersing_etal(u_ll, u_rr, orientation::Integer,
                                  equations::ShallowWaterExnerEquations2D)
    # Unpack left and right state
    _, h_v1_ll, h_v2_ll, _ = u_ll
    _, h_v1_rr, h_v2_rr, _ = u_rr

    # Get the velocities on either side
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Get the sediment discharge on either side
    q_s1_ll, q_s2_ll = sediment_discharge(u_ll, equations)
    q_s1_rr, q_s2_rr = sediment_discharge(u_rr, equations)

    # Average each factor of products in flux
    h_v1_avg = 0.5f0 * (h_v1_ll + h_v1_rr)
    h_v2_avg = 0.5f0 * (h_v2_ll + h_v2_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    q_s1_avg = 0.5f0 * (q_s1_ll + q_s1_rr)
    q_s2_avg = 0.5f0 * (q_s2_ll + q_s2_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = h_v1_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = q_s1_avg
    else # orientation == 2
        f1 = h_v2_avg
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = q_s2_avg
    end

    return SVector(f1, f2, f3, f4)
end

"""
    dissipation_roe(u_ll, u_rr, orientation::Integer,
                    equations::ShallowWaterExnerEquations2D)
Roe-type dissipation term for the [`ShallowWaterExnerEquations2D`](@ref) with an approximate Roe average
for the sediment discharge `q_s`.
"""
@inline function dissipation_roe(u_ll, u_rr, orientation::Integer,
                                 equations::ShallowWaterExnerEquations2D)
    r = equations.r
    g = equations.gravity
    z = zero(eltype(u_ll))

    # Get the velocities
    v1_ll, v2_ll = velocity(u_ll, equations)
    v1_rr, v2_rr = velocity(u_rr, equations)

    # Compute approximate Roe averages.
    # The actual Roe average for the sediment height `h_b` depends on the sediment and
    # friction model and an explicit formula is not always available.
    # Therefore we only use an approximation here.
    h_avg = 0.5f0 * (u_ll[1] + u_rr[1])
    v1_avg = (sqrt(u_ll[1]) * v1_ll + sqrt(u_rr[1]) * v1_rr) /
             (sqrt(u_ll[1]) + sqrt(u_rr[1]))
    v2_avg = (sqrt(u_ll[1]) * v2_ll + sqrt(u_rr[1]) * v2_rr) /
             (sqrt(u_ll[1]) + sqrt(u_rr[1]))
    h_b_avg = 0.5f0 * (u_ll[4] + u_rr[4])

    # Compute the nontrivial eigenvalues using Cardano's formula
    # The known eigenvalue of `v1` or `v2` associated with the contact wave is returned last.
    λ1, λ2, λ3, λ4 = eigvals_cardano(SVector(h_avg, h_avg * v1_avg, h_avg * v2_avg,
                                             h_b_avg),
                                     orientation, equations)

    # Compute the effective sediment height at the averaged solution state
    u_avg = SVector(h_avg, h_avg * v1_avg, h_avg * v2_avg, h_b_avg)
    h_s_avg = effective_sediment_height(u_avg, equations)

    # Build the right eigenvector matrix and its inverse in the appropriate direction
    if orientation == 1
        # Compute gradients of q_s1_avg using automatic differentiation.
        # Introduces a closure to make q_s2_avg a function of u only. This is necessary since the
        # gradient function only accepts functions of one variable.
        dq_s_dh, dq_s_dhv1, dq_s_dhv2, _ = Trixi.ForwardDiff.gradient(u -> sediment_discharge(u,
                                                                                              equations)[1],
                                                                      u_avg)

        # Precompute some common expressions
        c1 = g * (h_avg + h_s_avg)
        c2 = g * (h_avg + h_s_avg / r)

        # Eigenvector matrix
        r41 = ((v1_avg - λ1)^2 - c1) / c2
        r42 = ((v1_avg - λ2)^2 - c1) / c2
        r43 = ((v1_avg - λ3)^2 - c1) / c2

        # Workaround to avoid division by zero if `dq_s_dhv2` is close to zero.
        if abs(dq_s_dhv2) > eps(eltype(u_ll))
            r34 = -(dq_s_dh + λ4 * (dq_s_dhv1 + c1 / c2)) / dq_s_dhv2
            R = @SMatrix [[1 1 1 1]; [λ1 λ2 λ3 λ4]; [v2_avg v2_avg v2_avg r34];
                          [r41 r42 r43 -c1 / c2]]
        else # dq_s_dhv2 ≈ 0
            R = @SMatrix [[1 1 1 z]; [λ1 λ2 λ3 z]; [v2_avg v2_avg v2_avg 1];
                          [r41 r42 r43 z]]
        end

        # Inverse eigenvector matrix
        d1 = (λ1 - λ2) * (λ1 - λ3)
        d2 = (λ2 - λ1) * (λ2 - λ3)
        d3 = (λ3 - λ2) * (λ3 - λ1)

        # Workaround to avoid division by zero if `dq_s_dhv2` is close to zero.
        if abs(dq_s_dhv2) > eps(eltype(u_ll))
            D = r34 - v2_avg
            r_inv11 = (v2_avg * dq_s_dhv2 * (v1_avg - λ2) * (v1_avg - λ3) +
                       (v1_avg^2 - λ2 * λ3 - c1) *
                       (v2_avg * dq_s_dhv2 + dq_s_dh + (dq_s_dhv1 + c1 / c2) * v1_avg)) /
                      (D * d1 * dq_s_dhv2)
            r_inv21 = (v2_avg * dq_s_dhv2 * (v1_avg - λ1) * (v1_avg - λ3) +
                       (v1_avg^2 - λ1 * λ3 - c1) *
                       (v2_avg * dq_s_dhv2 + dq_s_dh + (dq_s_dhv1 + c1 / c2) * v1_avg)) /
                      (D * d2 * dq_s_dhv2)
            r_inv31 = (v2_avg * dq_s_dhv2 * (v1_avg - λ1) * (v1_avg - λ2) +
                       (v1_avg^2 - λ1 * λ2 - c1) *
                       (v2_avg * dq_s_dhv2 + dq_s_dh + (dq_s_dhv1 + c1 / c2) * v1_avg)) /
                      (D * d3 * dq_s_dhv2)
            R_inv = @SMatrix [r_inv11 (2 * v1_avg - λ2 - λ3)/d1 -(v1_avg - λ2) * (v1_avg - λ3)/(D * d1) c2/d1;
                              r_inv21 (2 * v1_avg - λ1 - λ3)/d2 -(v1_avg - λ1) * (v1_avg - λ3)/(D * d2) c2/d2;
                              r_inv31 (2 * v1_avg - λ2 - λ1)/d3 -(v1_avg - λ1) * (v1_avg - λ2)/(D * d3) c2/d3;
                              -v2_avg/D 0 1/D 0]
        else # dq_s_dhv2 ≈ 0
            R_inv = @SMatrix [(c1 - v1_avg^2 + λ2 * λ3)/d1 (2 * v1_avg - λ2 - λ3)/d1 z c2/d1;
                              (c1 - v1_avg^2 + λ1 * λ3)/d2 (2 * v1_avg - λ1 - λ3)/d2 z c2/d2;
                              (c1 - v1_avg^2 + λ1 * λ2)/d3 (2 * v1_avg - λ2 - λ1)/d3 z c2/d3;
                              -v2_avg z 1 z]
        end
    else # orientation == 2
        # Compute gradients of q_s2_avg using automatic differentiation.
        # Introduces a closure to make q_s2_avg a function of u only. This is necessary since the
        # gradient function only accepts functions of one variable.
        dq_s_dh, dq_s_dhv1, dq_s_dhv2, _ = Trixi.ForwardDiff.gradient(u -> sediment_discharge(u,
                                                                                              equations)[2],
                                                                      u_avg)

        # Precompute some common expressions
        c1 = g * (h_avg + h_s_avg)
        c2 = g * (h_avg + h_s_avg / r)

        # Eigenvector matrix
        r41 = ((v2_avg - λ1)^2 - c1) / c2
        r42 = ((v2_avg - λ2)^2 - c1) / c2
        r43 = ((v2_avg - λ3)^2 - c1) / c2

        # Workaround to avoid division by zero if `dq_s_dhv1` is close to zero.
        if abs(dq_s_dhv1) > eps(eltype(u_ll))
            r24 = -(dq_s_dh + λ4 * (dq_s_dhv2 + c1 / c2)) / dq_s_dhv1
            R = @SMatrix [[1 1 1 1]; [v1_avg v1_avg v1_avg r24]; [λ1 λ2 λ3 λ4];
                          [r41 r42 r43 -c1 / c2]]
        else # dq_s_dhv1 ≈ 0
            R = @SMatrix [[1 1 1 z]; [v1_avg v1_avg v1_avg 1]; [λ1 λ2 λ3 z];
                          [r41 r42 r43 z]]
        end

        # Inverse eigenvector matrix
        d1 = (λ1 - λ2) * (λ1 - λ3)
        d2 = (λ2 - λ1) * (λ2 - λ3)
        d3 = (λ3 - λ2) * (λ3 - λ1)

        # Workaround to avoid division by zero if `dq_s_dhv1` is close to zero.
        if abs(dq_s_dhv1) > eps(eltype(u_ll))
            D = r24 - v1_avg
            r_inv11 = (v1_avg * dq_s_dhv1 * (v2_avg - λ2) * (v2_avg - λ3) +
                       (v2_avg^2 - λ2 * λ3 - c1) *
                       (v1_avg * dq_s_dhv1 + dq_s_dh + (dq_s_dhv2 + c1 / c2) * v2_avg)) /
                      (D * d1 * dq_s_dhv1)
            r_inv21 = (v1_avg * dq_s_dhv1 * (v2_avg - λ1) * (v2_avg - λ3) +
                       (v2_avg^2 - λ1 * λ3 - c1) *
                       (v1_avg * dq_s_dhv1 + dq_s_dh + (dq_s_dhv2 + c1 / c2) * v2_avg)) /
                      (D * d2 * dq_s_dhv1)
            r_inv31 = (v1_avg * dq_s_dhv1 * (v2_avg - λ1) * (v2_avg - λ2) +
                       (v2_avg^2 - λ1 * λ2 - c1) *
                       (v1_avg * dq_s_dhv1 + dq_s_dh + (dq_s_dhv2 + c1 / c2) * v2_avg)) /
                      (D * d3 * dq_s_dhv1)
            R_inv = @SMatrix [r_inv11 -(v2_avg - λ2) * (v2_avg - λ3)/(D * d1) (2 * v2_avg - λ2 - λ3)/d1 c2/d1;
                              r_inv21 -(v2_avg - λ1) * (v2_avg - λ3)/(D * d2) (2 * v2_avg - λ1 - λ3)/d2 c2/d2;
                              r_inv31 -(v2_avg - λ1) * (v2_avg - λ2)/(D * d3) (2 * v2_avg - λ2 - λ1)/d3 c2/d3;
                              -v1_avg/D 1/D 0 0]
        else # dq_s_dhv1 ≈ 0
            R_inv = @SMatrix [(c1 - v2_avg^2 + λ2 * λ3)/d1 z (2 * v2_avg - λ2 - λ3)/d1 c2/d1;
                              (c1 - v2_avg^2 + λ1 * λ3)/d2 z (2 * v2_avg - λ1 - λ3)/d2 c2/d2;
                              (c1 - v2_avg^2 + λ1 * λ2)/d3 z (2 * v2_avg - λ2 - λ1)/d3 c2/d3;
                              -v1_avg 1 z z]
        end
    end

    # Eigenvalue absolute value matrix
    Λ_abs = @SMatrix [abs(λ1) z z z; z abs(λ2) z z; z z abs(λ3) z; z z z abs(λ4)]

    # Compute dissipation
    diss = SVector(-0.5f0 * R * Λ_abs * R_inv * (u_rr - u_ll))

    return SVector(diss[1], diss[2], diss[3], diss[4])
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterExnerEquations2D)
    return max(maximum(abs, eigvals_cardano(u_rr, orientation, equations)),
               maximum(abs, eigvals_cardano(u_ll, orientation, equations)))
end

@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterExnerEquations2D)
    lambdas_x = eigvals_cardano(u, 1, equations)
    lambdas_y = eigvals_cardano(u, 2, equations)
    return maximum(abs, lambdas_x), maximum(abs, lambdas_y)
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterExnerEquations2D)
    h, hv1, hv2, _ = u

    v1 = hv1 / h
    v2 = hv2 / h

    return SVector(v1, v2)
end

"""
    effective_sediment_height(u, equations::ShallowWaterExnerEquations2D)

Compute the "effective" sediment height `h_s` of the sediment discharge `q_s = h_s v_{1,2}`
for the Grass model.
Note, the inverse porosity scaling is put onto this quantity as a design decision.

- F. Benkhaldoun, S. Sahmim, M. Seaïd (2010)
  A two-dimensional finite volume morphodynamic model on unstructured triangular grids
  [DOI: 10.1002/fld.2129](https://doi.org/10.1002/fld.2129)
"""
@inline function effective_sediment_height(u,
                                           equations::ShallowWaterExnerEquations2D{T, S,
                                                                                   GrassModel{T}}) where {
                                                                                                          T,
                                                                                                          S
                                                                                                          }
    (; A_g, m_g) = equations.sediment_model
    v1, v2 = velocity(u, equations)
    v_norm = sqrt(v1^2 + v2^2)

    return equations.porosity_inv * A_g * v_norm^(m_g - 1)
end

"""
    effective_sediment_height(u, equations::ShallowWaterExnerEquations2D)

Compute the "effective" sediment height `h_s` of the sediment discharge `q_s = h_s v_{1,2}`
for Shields stress models.
This 2D extension is based upon the Meyer-Peter-Müller model described in the given reference.
Note, the inverse porosity scaling is put onto this quantity as a design decision.

- M.J. Castro Díaz, E.D. Fernández-Nieto, A.M. Ferreiro, and C. Parés (2009)\
  Two-dimensional sediment transport models in shallow water equations.
  A second order finite volume approach on unstructured meshes\
  [DOI: 10.1016/j.cma.2009.03.001](https://doi.org/10.1016/j.cma.2009.03.001)
"""
@inline function effective_sediment_height(u,
                                           equations::ShallowWaterExnerEquations2D{T, S,
                                                                                   ShieldsStressModel{T}}) where {
                                                                                                                  T,
                                                                                                                  S
                                                                                                                  }
    (; gravity, rho_f, rho_s, porosity_inv) = equations
    (; m_1, m_2, m_3, k_1, k_2, k_3, theta_c, d_s) = equations.sediment_model

    tau1, tau2 = shear_stress(u, equations)

    # Shields stress coinciding with the velocity vector
    theta = rho_f * sqrt(tau1^2 + tau2^2) / (gravity * (rho_s - rho_f) * d_s)

    # Characteristic discharge
    Q = d_s * sqrt(gravity * (rho_s / rho_f - 1) * d_s)

    # Compute the effective sediment height with workaround to avoid division by zero
    v1, v2 = velocity(u, equations)
    v_norm = sqrt(v1^2 + v2^2)
    h_s = zero(eltype(u))
    if v_norm > eps(eltype(u))
        h_s = (porosity_inv * Q / v_norm * k_1 * theta^m_1 *
               (max(theta - k_2 * theta_c, 0))^m_2 *
               (max(sqrt(theta) - k_3 * sqrt(theta_c), 0))^m_3)
    end

    return h_s
end

# Shear stress formulation using a coefficient to take into account different friction models
@inline function shear_stress(u, equations::ShallowWaterExnerEquations2D)
    v1, v2 = velocity(u, equations)
    v_norm = sqrt(v1^2 + v2^2)
    g = equations.gravity
    shear_coeff = shear_stress_coefficient(u, equations.friction)

    return SVector(g * shear_coeff * v1 * v_norm, g * shear_coeff * v2 * v_norm)
end

# TODO: how would this work in a normal direction?
#       handled inside the flux to compute the `q_s` in the normal direction
"""
    sediment_discharge(u, equations::ShallowWaterExnerEquations2D)

Compute the sediment discharge `q_s = h_s v_{1,2}` for a generic sediment model.
The dependency on the sediment model, like Grass or Shields, is inside [`effective_sediment_height`](@ref).
"""
@inline function sediment_discharge(u, equations::ShallowWaterExnerEquations2D)
    h_s = effective_sediment_height(u, equations)
    v1, v2 = velocity(u, equations)

    return SVector(h_s * v1, h_s * v2)
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterExnerEquations2D)
    h, _, _, h_b = u

    H = h + h_b
    v1, v2 = velocity(u, equations)

    return SVector(H, v1, v2, h_b)
end

# Convert conservative variables to entropy variables
@inline function Trixi.cons2entropy(u, equations::ShallowWaterExnerEquations2D)
    h, _, _, h_b = u
    (; gravity, r) = equations

    v1, v2 = velocity(u, equations)

    w1 = r * (gravity * (h + h_b) - 0.5f0 * (v1^2 + v2^2))
    w2 = r * v1
    w3 = r * v2
    w4 = gravity * (r * h + h_b)

    return SVector(w1, w2, w3, w4)
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::ShallowWaterExnerEquations2D)
    H, v1, v2, h_b = prim

    h = H - h_b
    hv1 = h * v1
    hv2 = h * v2

    return SVector(h, hv1, hv2, h_b)
end

@inline function Trixi.waterheight(u, equations::ShallowWaterExnerEquations2D)
    return u[1]
end

# Indicator variable used for the shock capturing in `IndicatorHennemannGassnerShallowWater`
@inline function water_sediment_height(u, equations::ShallowWaterExnerEquations2D)
    return equations.gravity * u[1] * u[4]
end

# Mathematical entropy function
@inline function Trixi.entropy(u, equations::ShallowWaterExnerEquations2D)
    h, _, _, h_b = u
    v1, v2 = velocity(u, equations)
    (; gravity, r) = equations

    return 0.5f0 * r * h * (v1^2 + v2^2) + 0.5f0 * gravity * (r * h^2 + h_b^2) +
           r * gravity * h * h_b
end

# Calculate the error for the "lake-at-rest" test case where H = h + h_b should
# be a constant value over time.
@inline function Trixi.lake_at_rest_error(u, equations::ShallowWaterExnerEquations2D)
    h, _, _, h_b = u
    return abs(equations.H0 - (h + h_b))
end

# Trigonometric version of Cardano's method to compute the nontrivial roots of a cubic polynomial
#   (x - v1,2)(x^3 + bx^2 + cx + d) = 0
# for the eigenvalues of the [`ShallowWaterExnerEquations2D`[(@ref)] flux Jacobian.
# This exploits that we know that either `v1` or `v2` is an eigenvalue (depending on the orientation)
# The eigenvalue that is equal to the velocity is associated with the contact wave
# in the Riemann fan and is returned as the last entry of the eigenvalue vector
# as expected by the `dissipation_roe`.
# TODO: Write specialized function for a similar strategy with `normal_direction` where
#       the characteristic polynomial is
#          (x - vn)(x^3 + bx^2 + cx + d) = 0
#       with vn = n1 v1 + n2 v2
# Note, assumes only real roots.
@inline function eigvals_cardano(u, orientation::Integer,
                                 equations::ShallowWaterExnerEquations2D)
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)
    g = equations.gravity
    r = equations.r

    # Compute the effective sediment height
    h_s = effective_sediment_height(u, equations)

    # Set the coefficients for the original cubic equation x^3 + bx^2 + cx + dx = 0
    # Note, some values change depending on orientation.
    if orientation == 1
        # Compute gradients of q_s1 using automatic differentiation.
        # Introduces a closure to make q_s1 a function of u only. This is necessary since the
        # gradient function only accepts functions of one variable.
        dq_s_dh, dq_s_dhv1, dq_s_dhv2, _ = Trixi.ForwardDiff.gradient(u -> sediment_discharge(u,
                                                                                              equations)[1],
                                                                      u)

        b = -2 * v1
        c = v1^2 - g * (h + h_s) - g * dq_s_dhv1 * (h + h_s / r)
        d = -g * dq_s_dh * (h + h_s / r) - g * v2 * dq_s_dhv2 * (h + h_s / r)
        # Set the known eigenvalue
        λ4 = v1
    else # orientation == 2
        # Compute gradients of q_s2 using automatic differentiation.
        # Introduces a closure to make q_s2 a function of u only. This is necessary since the
        # gradient function only accepts functions of one variable.
        dq_s_dh, dq_s_dhv1, dq_s_dhv2, _ = Trixi.ForwardDiff.gradient(u -> sediment_discharge(u,
                                                                                              equations)[2],
                                                                      u)

        b = -2 * v2
        c = v2^2 - g * (h + h_s) - g * dq_s_dhv2 * (h + h_s / r)
        d = -g * dq_s_dh * (h + h_s / r) - g * v1 * dq_s_dhv1 * (h + h_s / r)
        # Set the known eigenvalue
        λ4 = v2
    end

    # Once coefficients are computed we apply the trigonometric Cardano method

    # Create the coefficients of the depressed cubic equation t^3 + pt + q = 0
    p = c - b^2 / 3
    q = 2 * b^3 / 27 - b * c / 3 + d

    # Avoid round-off errors
    theta = clamp(3 * q / (2 * p) * sqrt(-3 / p), -1.0, 1.0)

    # Save common (but expensive) terms in the cubic root formula
    phi = acos(theta) / 3
    coeff = 2 * sqrt(-p / 3)
    shift = -b / 3

    # Use trigonometric form of Cardano to compute the three roots
    # Note, the fourth eigenvalue λ4 is set above in the if statement
    λ1 = shift + coeff * cos(phi)
    λ2 = shift + coeff * cos(phi - 2 * π / 3)
    λ3 = shift + coeff * cos(phi - 4 * π / 3)

    return SVector(λ1, λ2, λ3, λ4)
end
end # @muladd
