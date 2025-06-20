module BackwardDifferenceFormula

using StaticArrays

@inline tupletoSMatrix(A::NTuple{N, NTuple{N, T}}) where {N,T} = SMatrix{N,N,T}(A[i][j] for i in 1:N, j in 1:N)

@generated function Vandermonde(dt::Union{SVector{order},NTuple{order}}) where {order}
    quote
        @inline 
        A = Base.Cartesian.@ntuple $order i -> begin
            _fact = 1 / factorial(i - 1)
            n     = i - 1
            Base.Cartesian.@ntuple $order j -> begin
                @fastmath dt[j]^n * _fact 
            end
        end
        tupletoSMatrix(A)
    end
end

"""
    coeff = bdf_coefficients(tshift)

Computes Backward Difference Formula coefficient for variable step sizes using a list of times values.
The reference current time should be 0, and the values of relative time for the `n` preceding steps should be provided.  

# Input
- `tshift`: reference time including previous steps.
     e.g., [-Δt, 0] generates BDF1 coefficients while  [-Δt0-Δt, -Δt, 0] produces BDF2 coefficients
# Output
    `coeff`` contains BDF coefficients. The first entry correspont to the next field value, e.g.:
    for BDF1:  (coeff[1]*u + coeff[2]*u0) = (u - u0) / Δt ≈ ∂u∂t
    for BDF2:  (coeff[1]*u + coeff[2]*u0 + coeff[3]*u00) ≈ ∂u∂t
"""
function bdf_coefficients(tshift::Union{SVector{order, T}, NTuple{order, T}}) where {order, T}
    # Construct Vandermonde matrix
    A = Vandermonde(tshift)

    # # Righ hand side
    rhs    = @MVector zeros(T, order)
    rhs[2] = one(T)  # Representing the derivative

    # Solve
    coeffs = A \ SVector(rhs)

    return reverse(coeffs)
end

export bdf_coefficients

end # module BackwardDifferenceFormula




