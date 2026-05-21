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
\partial_t h_b + \partial_x q_s + \partial_y q_s = 0,
\end{cases}
```
The unknown quantities are the water and sediment height ``h``, ``h_b`` and the velocities ``v_1``, ``v_2``.
The sediment discharges ``q_s1(h, hv)`` and ``q_s1(h, hv)`` are determined by the `sediment_model` and is used to determine
the active sediment heights ``h_s1 = q_s1 / v_1`` and ``h_s2 = q_s2 / v_2``.
Furthermore ``\tau`` denotes the shear stress at the water-sediment interface and is determined by
the `friction` model.
The gravitational acceleration is denoted by ``g``, and ``\rho_f`` and ``\rho_s`` are the fluid and sediment
densities, respectively. The density ratio is given by ``r = \rho_f / \rho_s``, where ``r`` lies between
``0 < r < 1`` as the fluid density ``\rho_f`` should be smaller than the sediment density ``\rho_s``.

The conservative variable water height ``h`` is measured from the sediment height ``h_b``, therefore
one also defines the total water height as ``H = h + h_b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water
height that is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

The entropy conservative formulation has been derived in the paper:
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

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_convergence_test(x, t, equations::ShallowWaterExnerEquations2D)

A smooth initial condition used for convergence tests in combination with
[`Trixi.source_terms_convergence_test`](@ref).
"""
@inline function Trixi.initial_condition_convergence_test(x, t,
                                                          equations::ShallowWaterExnerEquations2D)
    ω = sqrt(2) * pi

    h = 2 + cos(ω * x[1]) * cos(ω * x[2]) * cos(ω * t)
    v1 = 0.5f0
    v2 = -0.65f0
    h_b = 2 + sin(ω * x[1]) * sin(ω * x[2]) * cos(ω * t)

    return SVector(h, h * v1, h * v2, h_b)
end

# """
#     source_terms_convergence_test(u, x, t, equations::ShallowWaterExnerEquations1D{T, S, GrassModel{T}}) where {T, S}

# Source terms used for convergence tests in combination with [`Trixi.initial_condition_convergence_test`](@ref)
# when using the the [`GrassModel`](@ref) model.

# To use this source term the equations must be set to:
# ```julia
# equations = ShallowWaterExnerEquations2D(gravity = 10.0, rho_f = 0.5,
#                                             rho_s = 1.0, porosity = 0.5,
#                                             friction = ManningFriction(n = 0.0),
#                                             sediment_model = GrassModel(A_g = 0.01))
# ```
# """
# # TODO: needs updated
# @inline function Trixi.source_terms_convergence_test(u, x, t,
#                                                      equations::ShallowWaterExnerEquations2D{T,
#                                                                                              S,
#                                                                                              GrassModel{T}}) where {
#                                                                                                                     T,
#                                                                                                                     S
#                                                                                                                     }
#     ω = sqrt(2) * pi
#     A_g = equations.sediment_model.A_g

#     h = -cos(x[1] * ω) * sin(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω
#     hv = -0.5f0 * cos(x[1] * ω) * sin(t * ω) * ω -
#          0.25f0 * sin(x[1] * ω) * cos(t * ω) * ω +
#          10 * A_g *
#          (cos(x[1] * ω) * cos(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω) +
#          10 * (2 + cos(x[1] * ω) * cos(t * ω)) *
#          (cos(x[1] * ω) * cos(t * ω) * ω - sin(x[1] * ω) * cos(t * ω) * ω)
#     h_b = -sin(x[1] * ω) * sin(t * ω) * ω
#     return SVector(h, hv, h_b)
# end

# """
#     source_terms_convergence_test(u, x, t, equations::ShallowWaterExnerEquations2D{T, S, ShieldsStressModel{T}}) where {T, S}

# Source terms used for convergence tests in combination with [`Trixi.initial_condition_convergence_test`](@ref)
# when using the [`MeyerPeterMueller`](@ref) model.

# To use this source term the equations must be set to:
# ```julia
# equations = ShallowWaterExnerEquations2D(gravity = 10.0, rho_f = 0.5,
#                                          rho_s = 1.0, porosity = 0.5,
#                                          friction = ManningFriction(n = 0.01),
#                                          sediment_model = MeyerPeterMueller(theta_c = 0.0,
#                                                                             d_s = 1e-3))
# ```
# """
# #TODO: needs updated
# @inline function Trixi.source_terms_convergence_test(u, x, t,
#                                                      equations::ShallowWaterExnerEquations2D{T,
#                                                                                              S,
#                                                                                              ShieldsStressModel{T}}) where {
#                                                                                                                             T,
#                                                                                                                             S
#                                                                                                                             }
#     ω = sqrt(2) * pi
#     (; gravity, porosity_inv, rho_f, rho_s, r) = equations

#     n = equations.friction.n

#     # Constant expression from the MPM model
#     c = sqrt(gravity * (1 / r - 1)) * 8 * porosity_inv *
#         (rho_f / (rho_s - rho_f))^(3 / 2) * n^3

#     h = -cos(x[1] * ω) * sin(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω

#     hv = ((5 * c *
#            (cos(x[1] * ω) * cos(t * ω) * ω - 0.5f0 * sin(x[1] * ω) * cos(t * ω) * ω)) /
#           ((2 + cos(x[1] * ω) * cos(t * ω))^0.5) -
#           0.5f0 * cos(x[1] * ω) * sin(t * ω) * ω -
#           0.25f0 * sin(x[1] * ω) * cos(t * ω) * ω +
#           10 * (2 + cos(x[1] * ω) * cos(t * ω)) *
#           (cos(x[1] * ω) * cos(t * ω) * ω - sin(x[1] * ω) * cos(t * ω) * ω))

#     h_b = ((0.5f0 * ((0.125f0 * c) / (2 + cos(x[1] * ω) * cos(t * ω))) * sin(x[1] * ω) *
#             cos(t * ω) * ω) / ((2 + cos(x[1] * ω) * cos(t * ω))^0.5) -
#            sin(x[1] * ω) * sin(t * ω) * ω)

#     return SVector(h, hv, h_b)
# end

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
    qs1, qs2 = q_s(u, equations)

    if orientation == 1
        f1 = hv1
        f2 = hv1 * v1
        f3 = hv2 * v1
        f4 = qs1
    else # orientation == 2
        f1 = hv2
        f2 = hv1 * v2
        f3 = hv2 * v2
        f4 = qs2
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

    # Calculate velocities
    v1_ll, v2_ll = velocity(u_ll, equations)

    # Get the local sediment discharges
    qs1_ll, qs2_ll = q_s(u_ll, equations)

    # Workaround to avoid division by zero, when computing the effective sediment height
    if orientation == 1
        if abs(v1_ll) < eps(typeof(h_ll))
            h_s_ll = 0
        else
            h_s_ll = qs1_ll / v1_ll
        end
    else # orientation == 2
        if abs(v2_ll) < eps(typeof(h_ll))
            h_s_ll = 0
        else
            h_s_ll = qs2_ll / v2_ll
        end
    end

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
                                     equations::ShallowWaterMultiLayerEquations2D)

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
    qs1_ll, qs2_ll = q_s(u_ll, equations)
    qs1_rr, qs2_rr = q_s(u_rr, equations)

    # Average each factor of products in flux
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)

    # Calculate fluxes depending on orientation
    if orientation == 1
        f1 = 0.5f0 * (h_v1_ll + h_v1_rr)
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = 0.5f0 * (qs1_ll + qs1_rr)
    else # orientation == 2
        f1 = 0.5f0 * (h_v2_ll + h_v2_rr)
        f2 = f1 * v1_avg
        f3 = f1 * v2_avg
        f4 = 0.5f0 * (qs2_ll + qs2_rr)
    end

    return SVector(f1, f2, f3, f4)
end

"""
    dissipation_roe(u_ll, u_rr, orientation,
                    equations::ShallowWaterExnerEquations2D)
Roe-type dissipation term for the [`ShallowWaterExnerEquations1D`](@ref) with an approximate Roe average
for the sediment discharge `q_s`.
"""
@inline function dissipation_roe(u_ll, u_rr, orientation,
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

    # Compute the eigenvalues using Ferrari's formula
    lambdas = eigvals_ferrari(SVector(h_avg, h_avg * v1_avg, h_avg * v2_avg, h_b_avg),
                              orientation, equations)

    # Compute the sediment discharge at the averaged state
    q_s1_avg, q_s2_avg = q_s(SVector(h_avg, h_avg * v1_avg, h_avg * v2_avg, h_b_avg),
                             equations)

    # Build the right eigenvector matrix and its inverse in the appropriate direction
    if orientation == 1
        # Workaround to avoid division by zero, when computing the effective sediment height
        if abs(v1_avg) < eps(typeof(h_avg))
            h_s_avg = z
        else
            h_s_avg = q_s1_avg / v1_avg
        end

        # Precompute some common expressions
        c1 = g * (h_avg + h_s_avg)
        c2 = g * (h_avg + h_s_avg / r)

        # Identify which eigenvalue is closest to `v1_avg` as it is associated
        # with a contact wave with a degenerate eigenvector
        tol = 1e-10 * max(1, abs(v1_avg), sqrt(c1))
        contact_idx = 1
        min_err = abs(lambdas[1] - v1_avg)
        for i in 2:4
            err = abs(lambdas[i] - v1_avg)
            if err < min_err
                min_err = err
                contact_idx = i
            end
        end

        # Swap the contact eigenvalue into the last slot
        lambdas_swapped = swap_contact_ev_to_last(lambdas, contact_idx)

        # Unpack the eigenvalues for convenience
        λ1 = lambdas_swapped[1]
        λ2 = lambdas_swapped[2]
        λ3 = lambdas_swapped[3]
        λ4 = lambdas_swapped[4]

        # Eigenvector matrix
        r41 = ((v1_avg - λ1)^2 - c1) / c2
        r42 = ((v1_avg - λ2)^2 - c1) / c2
        r43 = ((v1_avg - λ3)^2 - c1) / c2
        R = @SMatrix [[1 1 1 z]; [λ1 λ2 λ3 z]; [v2_avg v2_avg v2_avg 1];
                      [r41 r42 r43 z]]

        # Inverse eigenvector matrix
        d1 = (λ1 - λ2) * (λ1 - λ3)
        d2 = (λ2 - λ1) * (λ2 - λ3)
        d3 = (λ3 - λ2) * (λ3 - λ1)
        R_inv = @SMatrix [(c1 - v1_avg^2 + λ2 * λ3)/d1 (2 * v1_avg - λ2 - λ3)/d1 z c2/d1;
                          (c1 - v1_avg^2 + λ1 * λ3)/d2 (2 * v1_avg - λ1 - λ3)/d2 z c2/d2;
                          (c1 - v1_avg^2 + λ1 * λ2)/d3 (2 * v1_avg - λ2 - λ1)/d3 z c2/d3;
                          -v2_avg z 1 z]
    else # orientation == 2
        # Workaround to avoid division by zero, when computing the effective sediment height
        if abs(v2_avg) < eps(typeof(h_avg))
            h_s_avg = z
        else
            h_s_avg = q_s2_avg / v2_avg
        end

        # Precompute some common expressions
        c1 = g * (h_avg + h_s_avg)
        c2 = g * (h_avg + h_s_avg / r)

        # Identify which eigenvalue is closest to `v2_avg` as it is associated
        # with a contact wave with a degenerate eigenvector
        tol = 1e-10 * max(1, abs(v2_avg), sqrt(c1))
        contact_idx = 1
        min_err = abs(lambdas[1] - v2_avg)
        for i in 2:4
            err = abs(lambdas[i] - v2_avg)
            if err < min_err
                min_err = err
                contact_idx = i
            end
        end

        # Swap the contact eigenvalue into the last slot
        lambdas_swapped = swap_contact_ev_to_last(lambdas, contact_idx)

        # Unpack the eigenvalues for convenience
        λ1 = lambdas_swapped[1]
        λ2 = lambdas_swapped[2]
        λ3 = lambdas_swapped[3]
        λ4 = lambdas_swapped[4]

        # Eigenvector matrix
        r41 = ((v2_avg - λ1)^2 - c1) / c2
        r42 = ((v2_avg - λ2)^2 - c1) / c2
        r43 = ((v2_avg - λ3)^2 - c1) / c2
        R = @SMatrix [[1 1 1 z]; [v1_avg v1_avg v1_avg 1]; [λ1 λ2 λ3 z];
                      [r41 r42 r43 z]]

        # Inverse eigenvector matrix
        d1 = (λ1 - λ2) * (λ1 - λ3)
        d2 = (λ2 - λ1) * (λ2 - λ3)
        d3 = (λ3 - λ2) * (λ3 - λ1)
        R_inv = @SMatrix [(c1 - v2_avg^2 + λ2 * λ3)/d1 z (2 * v2_avg - λ2 - λ3)/d1 c2/d1;
                          (c1 - v2_avg^2 + λ1 * λ3)/d2 z (2 * v2_avg - λ1 - λ3)/d2 c2/d2;
                          (c1 - v2_avg^2 + λ1 * λ2)/d3 z (2 * v2_avg - λ2 - λ1)/d3 c2/d3;
                          -v1_avg 1 z z]
    end

    # Eigenvalue absolute value matrix
    Λ_abs = @SMatrix [abs(λ1) z z z; z abs(λ2) z z; z z abs(λ3) z; z z z abs(λ4)]

    # Compute dissipation
    diss = SVector(-0.5f0 * R * Λ_abs * R_inv * (u_rr - u_ll))

    return SVector(diss[1], diss[2], diss[3], diss[4])
end

@inline function swap_contact_ev_to_last(lambdas, contact_idx)
    if contact_idx == 4
        return lambdas
    elseif contact_idx == 1
        return SVector(lambdas[2], lambdas[3], lambdas[4], lambdas[1])
    elseif contact_idx == 2
        return SVector(lambdas[1], lambdas[3], lambdas[4], lambdas[2])
    else # contact_idx == 3
        return SVector(lambdas[1], lambdas[2], lambdas[4], lambdas[3])
    end
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                           equations::ShallowWaterExnerEquations2D)
    return max(maximum(abs, eigvals_ferrari(u_rr, orientation, equations)),
               maximum(abs, eigvals_ferrari(u_ll, orientation, equations)))
end

@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterExnerEquations2D)
    lambdas_x = eigvals_ferrari(u, 1, equations)
    lambdas_y = eigvals_ferrari(u, 2, equations)
    return maximum(abs, lambdas_x), maximum(abs, lambdas_y)
end

#Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterExnerEquations2D)
    h, hv1, hv2, _ = u

    v1 = hv1 / h
    v2 = hv2 / h

    return SVector(v1, v2)
end

# Compute the sediment discharge for Shields stress models
# TODO: double check (somewhere?) that this 2d generalization makes sense
# TODO: how would this work in a normal direction?
@inline function q_s(u,
                     equations::ShallowWaterExnerEquations2D{T, S,
                                                             ShieldsStressModel{T}}) where {
                                                                                            T,
                                                                                            S
                                                                                            }
    (; gravity, rho_f, rho_s, porosity_inv) = equations
    (; m_1, m_2, m_3, k_1, k_2, k_3, theta_c, d_s) = equations.sediment_model

    tau1, tau2 = shear_stress(u, equations)
    theta1 = rho_f * abs(tau1) / (gravity * (rho_s - rho_f) * d_s)  # Shields stress in x-direction
    theta2 = rho_f * abs(tau2) / (gravity * (rho_s - rho_f) * d_s)  # Shields stress in y-direction

    Q = d_s * sqrt(gravity * (rho_s / rho_f - 1) * d_s) # Characteristic discharge

    return SVector((porosity_inv * Q * sign(theta1) * k_1 * theta1^m_1 *
                    (max(theta1 - k_2 * theta_c, 0))^m_2 *
                    (max(sqrt(theta1) - k_3 * sqrt(theta_c), 0))^m_3),
                   (porosity_inv * Q * sign(theta2) * k_1 * theta2^m_1 *
                    (max(theta2 - k_2 * theta_c, 0))^m_2 *
                    (max(sqrt(theta2) - k_3 * sqrt(theta_c), 0))^m_3))
end

# Compute the sediment discharge for the Grass model
@inline function q_s(u,
                     equations::ShallowWaterExnerEquations2D{T, S, GrassModel{T}}) where {
                                                                                          T,
                                                                                          S
                                                                                          }
    (; porosity_inv, sediment_model) = equations
    v1, v2 = velocity(u, equations)
    v_norm = sqrt(v1^2 + v2^2)
    return SVector(porosity_inv * sediment_model.A_g * v1 *
                   v_norm^(sediment_model.m_g - 1),
                   porosity_inv * sediment_model.A_g * v2 *
                   v_norm^(sediment_model.m_g - 1))
end

# Shear stress formulation using a coefficient to take into account different friction models
@inline function shear_stress(u, equations::ShallowWaterExnerEquations2D)
    v1, v2 = velocity(u, equations)
    v_norm = sqrt(v1^2 + v2^2)
    g = equations.gravity
    shear_coeff = shear_stress_coefficient(u, equations.friction)
    return SVector(g * shear_coeff * v1 * v_norm, g * shear_coeff * v2 * v_norm)
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

# Branch aware modified Ferrari formula to compute the roots of a quartic polynomial
#   x^4 + bx^3 + cx^2 + dx + e = 0
# in order to compute the eigenvalues of the [`ShallowWaterExnerEquations2D`[(@ref)].
# Broad outline of the procedure below:
#   1) Depress the quartic
#   2) Use Trigonometric version of Cardano's method roots of a resolvant cubic polynomial
#   3) Compute the stable square root
#   4) Compute the (paired) roots
# Note, assumes only real roots.
@inline function eigvals_ferrari(u, orientation::Integer,
                                 equations::ShallowWaterExnerEquations2D)
    h = waterheight(u, equations)
    v1, v2 = velocity(u, equations)
    g = equations.gravity
    r = equations.r

    if orientation == 1
        # Workaround to avoid division by zero, when computing the effective sediment height
        if abs(v1) > eps(typeof(h))
            q_s1, _ = q_s(u, equations)
            h_s1 = q_s1 / v1
            # Compute gradients of q_s using automatic differentiation.
            # Introduces a closure to make q_s a function of u only. This is necessary since the
            # gradient function only accepts functions of one variable.
            dq_s1_dh, dq_s1_dhv1, dq_s1_dhv2, _ = Trixi.ForwardDiff.gradient(u -> q_s1,
                                                                             u)
        else
            h_s1 = 0
            dq_s1_dh = 0
            dq_s1_dhv1 = 0
            dq_s1_dhv2 = 0
        end
        # Coefficients for the original quartic equation x^4 + bx^3 + cx^2 + dx + e
        b = -3 * v1
        c = 3 * v1^2 - g * (h + h_s1 / r) * dq_s1_dhv1 - g * (h + h_s1)
        d = -v1^3 + g * (h + h_s1) * v1 +
            g * (h + h_s1 / r) * (-dq_s1_dh + dq_s1_dhv1 * v1 - dq_s1_dhv2 * v2)
        e = g * (h + h_s1 / r) * v1 * (dq_s1_dh + v2 * dq_s1_dhv2)
    else # orientation == 2
        # Workaround to avoid division by zero, when computing the effective sediment height
        if abs(v2) > eps(typeof(h))
            _, q_s2 = q_s(u, equations)
            h_s2 = q_s2 / v2
            # Compute gradients of q_s using automatic differentiation.
            # Introduces a closure to make q_s a function of u only. This is necessary since the
            # gradient function only accepts functions of one variable.
            dq_s2_dh, dq_s2_dhv1, dq_s2_dhv2, _ = Trixi.ForwardDiff.gradient(u -> q_s2,
                                                                             u)
        else
            h_s2 = 0
            dq_s2_dh = 0
            dq_s2_dhv1 = 0
            dq_s2_dhv2 = 0
        end
        b = -3 * v2
        c = 3 * v2^2 - g * (h + h_s2 / r) * dq_s2_dhv2 - g * (h + h_s2)
        d = -v2^3 + g * (h + h_s2) * v2 +
            g * (h + h_s2 / r) * (-dq_s2_dh - dq_s2_dhv1 * v1 + dq_s2_dhv2 * v2)
        e = g * (h + h_s2 / r) * v2 * (dq_s2_dh + v1 * dq_s2_dhv1)
    end

    # Once coefficients are computed we apply the modified Ferrari method

    # Coefficients of the depressed quartic equation y^4 + py^2 + qy + r = 0
    # where λ = y - b/4
    p = c - 3 * b^2 / 8
    q = b^3 / 8 - b * c / 2 + d
    r = -3 * b^4 / 256 + b^2 * c / 16 - b * d / 4 + e

    # Next, we need the largest root of the resolvant cubic equation
    #  8m^3 - 4pm^2 - 8rm + 4pr - q^2 = 0
    # which is computed from the trigonometric version of Cardano's formula
    B = -p / 2
    C = -r
    D = (4 * p * r - q^2) / 8

    # Coefficient of the depressed cubic equation t^3 + Pt + Q = 0
    P = C - B^2 / 3
    Q = 2 * B^3 / 27 - B * C / 3 + D

    # Avoid round-off errors
    theta = clamp(3 * Q / (2 * P) * sqrt(-3 / P), -1.0, 1.0)

    # Save common (but expensive) terms in the cubic roots
    phi = acos(theta) / 3
    coeff = 2 * sqrt(-P / 3)

    # Use trigonometric form of Cardano to compute the three roots
    λ1_c = -B / 3 + coeff * cos(phi)
    λ2_c = -B / 3 + coeff * cos(phi - 2 * π / 3)
    λ3_c = -B / 3 + coeff * cos(phi - 4 * π / 3)

    # Take the maximum root of the resolvant equation
    m = max(λ1_c, λ2_c, λ3_c)

    # Compute stable square root (avoids cancellation in nested sqrt)
    R = sqrt(max(0, 2 * m - p))
    t = R < eps(typeof(h)) ? 0 : 2 * q / R
    D1 = sqrt(max(0, -(2 * m + p) - t))
    D2 = sqrt(max(0, -(2 * m + p) + t))

    # Roots of the original cubic equation
    λ1 = -b / 4 + 0.5f0 * (R + D1)
    λ2 = -b / 4 + 0.5f0 * (R - D1)
    λ3 = -b / 4 + 0.5f0 * (-R + D2)
    λ4 = -b / 4 + 0.5f0 * (-R - D2)

    return SVector(λ1, λ2, λ3, λ4)
end
end # @muladd
