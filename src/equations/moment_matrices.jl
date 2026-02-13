# Create the B tensor for the moment equations in a (somewhat) clever way.
# We use recursion relations for the Legendre polynomials and their derivative
# to build the necessary components of the integral
#
#   B_ijk = (2i + 1) integral_0^1 phi'_i (integral_0^zeta phi_j ds) phi_k  d zeta
#   i, j, k = 1, ..., N
#
# where phi'_i is a polynomial of degree at most `N-1`
# and \phi_j and \phi_k are polynomials of degree at most `N`.
# Therefore, `integral_0^zeta phi_j ds` is a polynomial of at most `N+1`.
# So, the highest degree of the integrand is a polynomial
# of degree `(N-1) + (N+1) + N = 3N`. We approximate the integrals with
# Legendre-Gauss quadrature which is exact for polynomials up to degree `2M+1`.
# Thus, if we select `M = ceil( (3N - 1) / 2 ) + 1` nodes, then the quadrature
# will be exact for the construction of the B tensor.
function compute_B_tensor(n::Int, RealT = Float64)
    # save preliminary values in pre-allocated memory
    B = Array{RealT, 3}(zeros(RealT, n, n, n))

    # Given the number of moments `n` determine how many Legendre-Gauss
    # nodes are required and compute them
    M = Int(ceil(1.5f0 * n - 0.5f0) + 1)
    lg_nodes, lg_weights = Trixi.gauss_nodes_weights(M, RealT)

    # Precompute the values of \phi and \phi' evaluated at the mapped LG nodes
    phi_values = Array{RealT, 2}(undef, n, M)
    phi_prime_values = Array{RealT, 2}(undef, n, M)
    for j in 1:n, m in 1:M
        x = 0.5f0 * (lg_nodes[m] + 1)
        phi_values[j, m], phi_prime_values[j, m] = shifted_legendre_polynomial_and_derivative(j,
                                                                                              x)
    end # j, m

    # Precompute all the internal integrals
    phi_j_integrals = Array{RealT, 2}(undef, n, M)
    compute_inner_integrals!(phi_j_integrals, lg_nodes, lg_weights)

    # Compute everything using the precomputed work arrays
    for i in 1:n, j in 1:n, k in 1:n
        # Legendre-Gauss quadrature loop for the outer integral
        res = 0
        for m in 1:M
            # Outer LG quadrature node and weight
            w_m = lg_weights[m]

            res += phi_prime_values[i, m] * phi_j_integrals[j, m] * phi_values[k, m] * w_m
        end # m
        # Scale by (2i+1) / 2 (extra 1/2 from mapping integral from [0, 1] -> [-1, 1])
        B[i, j, k] = 0.5f0 * (2 * i + 1) * res
    end # i, j, k

    return Array{RealT, 3}(B)
end

# Directly compute the A tensor from B using the relation
#  A_ijk / (2i + 1) = - B_ijk / (2i + 1) - B_kji / (2k + 1)
function compute_A_tensor(B::Array{RealT, 3}) where {RealT}
    N = size(B, 1) # Get number of moments
    A = Array{RealT, 3}(zeros(RealT, N, N, N))  # initialize the A tensor

    for i in 1:N, j in 1:N, k in 1:N
        A[i, j, k] = -(2i + 1) * (B[i, j, k] / (2i + 1) + B[k, j, i] / (2k + 1))
    end # i, j, k

    return A
end

# Helper function to precompute and store every instance of the interior integral
# when computing the B tensor
#
#   integral_0^zeta phi_j(s) ds
#       = integral_0^(0.5*(xi_m + 1)) phi_j(s) d s
#       = 0.25 * (xi_m + 1) integral_-1^1 phi_j(0.25 (xi_m + 1) (eta_p + 1)) d eta
#
function compute_inner_integrals!(phi_j_integrals, lg_nodes, lg_weights)
    # save preliminary values in pre-allocated memory
    fill!(phi_j_integrals, 0)

    n, M = size(phi_j_integrals)

    for j in 1:n, m in 1:M
        # Outer LG quadrature node and weight
        xi_m = lg_nodes[m]

        # For a given `m` and `j` first compute the inner integral
        # Technically, this integral could use fewer LG nodes, but ignore for now
        for p in 1:M
            # Inner LG quadrature node and weight
            eta_p = lg_nodes[p]
            w_p = lg_weights[p]

            # Compute \phi_j at the (mapped) point 0.25 (xi_m + 1) (eta_p + 1)
            x = 0.25f0 * (xi_m + 1) * (eta_p + 1)
            phi_j_value, _ = shifted_legendre_polynomial_and_derivative(j, x)

            # Add result into the sum
            phi_j_integrals[j, m] += phi_j_value * w_p
        end # p
        # Jacobian from integration moving form [0, (xi_m+1)/2] -> [-1, 1]
        phi_j_integrals[j, m] *= 0.25f0 * (xi_m + 1)
    end # j, m

    return nothing
end

# Create the C matrix for the moment equations if one wants to include friction.
# We use recursion relations for the Legendre polynomial derivatives
# to build the necessary components of the integral
#
#   C_ij = (2i + 1) integral_0^1 phi'_i phi'_j  d zeta
#   i, j = 1, ..., N
#
# where phi'_i and \phi'_j are polynomials of degree at most `N-1`.
# So, the highest degree of the integrand is a polynomial
# of degree `(N-1) + (N-1) = 2N - 2`. We approximate the integrals with
# Legendre-Gauss quadrature which is exact for polynomials up to degree `2M+1`.
# Thus, if we select `M = ceil( (2N - 3) / 2 ) + 1` nodes, then the quadrature
# will be exact for the construction of the C matrix.
#
# It is worth noting that the shifted and unnormalized Legendre polynomials
# are no longer orthogonal on [0, 1], so this matrix is not diagonal.
# If one uses a normalized basis then we get orthogonality back as is a classic
# result for the Legendre polynomials also on [-1, 1], e.g.,
# Abramowitz & Stegun, Handbook of Mathematical Functions, Chapter 22
function compute_C_matrix(n::Int, RealT = Float64)
    # save preliminary values in pre-allocated memory
    C = Array{RealT, 2}(zeros(RealT, n, n))

    # Given the number of moments `n` determine how many Legendre-Gauss
    # nodes are required and compute them
    M = Int(ceil(n - 1.5f0) + 1)
    lg_nodes, lg_weights = Trixi.gauss_nodes_weights(M, RealT)

    # Compute everything on the fly (might be slow)
    for i in 1:n, j in 1:n
        # Legendre-Gauss quadrature loop for the outer integral
        res = 0
        for m in 1:M
            # Outer LG quadrature node and weight
            xi_m = lg_nodes[m]
            w_m = lg_weights[m]

            # Compute phi'_i and \phi'_j at the mapped node 0.5 (xi_m + 1)
            x = 0.5f0 * (xi_m + 1)
            _, phi_i_prime = shifted_legendre_polynomial_and_derivative(i, x)
            _, phi_j_prime = shifted_legendre_polynomial_and_derivative(j, x)

            # Add result into the sum
            res += phi_i_prime * phi_j_prime * w_m
        end # m
        # Scale by (2i+1) / 2 (extra 1/2 from mapping integral from [0, 1] -> [-1, 1])
        C[i, j] = 0.5f0 * (2 * i + 1) * res
    end # i, j

    return C
end

# The tensor B and matrix C are built from the shifted and unnormalized Legendre polynomials.
# They can be created from the Rodrigues' formula
#
#   phi_j = 1/j! d^j / dx^j (x - x^2)^j
#
# for j = 1, ..., N and x in [0, 1]. A more convenient way to build the polynomials
# is from the three term recurrence relation
#
#   phi_0 = 1,
#   phi_1 = 1 - 2x,
#   phi_j = ((2j-1) (1 - 2x) phi_{j-1} - (j-1) phi_{j-2}) / (j), j = 2, ..., N
#
# The derivatives of the shifted and unnormalized Legendre polynomials can also
# be built from a three term recurrence relation
#
#   phi'_0 = 0,
#   phi'_1 = -2,
#   phi'_j = phi'_{j-2} - 2 (2j-1) phi_{j-1}, j = 2, ..., N
#
# This routine computes and returns phi_N and phi'_N evaluated at a point `x`.
function shifted_legendre_polynomial_and_derivative(N::Int, x::Real)
    if N == 0
        poly = 1
        deriv = 0
    elseif N == 1
        poly = 1 - 2 * x
        deriv = -2
    else
        poly_Nm2 = 1
        poly_Nm1 = 1 - 2 * x
        deriv_Nm2 = 0
        deriv_Nm1 = -2

        poly = 0
        deriv = 0
        for j in 2:N
            poly = ((2 * j - 1) * (1 - 2 * x) * poly_Nm1 - (j - 1) * poly_Nm2) / j
            deriv = deriv_Nm2 - 2 * (2 * j - 1) * poly_Nm1
            poly_Nm2 = poly_Nm1
            poly_Nm1 = poly
            deriv_Nm2 = deriv_Nm1
            deriv_Nm1 = deriv
        end
    end

    return poly, deriv
end
