"""

    addtodiag!(A, v)

Efficiently add the vector `v` in-place to the main diagonal of the matrix `A` (overwriting
the latter). `length(v)` must match one of the dimensions of `A`, with the other dimension
being at least as large, and the element types of `v` must be promotable to those of `A`.
If so, yield the same result as `A .= diagm(v) .+ A`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> addtodiag!(A, v)
3×3 Array{Float64,2}:
 2.0  2.0   3.0
 2.0  6.0   6.0
 3.0  6.0  12.0
```

See also: [`addtodiag`](@ref), [`subtractfromdiag!`](@ref), [`subtractfromdiag`](@ref)

"""
function addtodiag!(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1})

    if (length(v) == size(A, 1) && length(v) <= size(A, 2))
        @inbounds for jj=1:size(A, 1)
            A[jj, jj] += v[jj];
        end
    elseif (length(v) == size(A, 2) && length(v) <= size(A, 1))
        @inbounds for jj=1:size(A, 2)
            A[jj, jj] += v[jj];
        end
    else
        error(string("addtodiag!(): Length of v (currently ", length(v),
                     ") must match at least one dimension of A (currently", size(A), ")"));
    end

end

"""

    addtodiag(A, v)

Efficiently add the vector `v` not in-place to the main diagonal of the matrix `A`,
returning the result. `length(v)` must match one of the dimensions of `A`, with the other
dimension being at least as large, and the element types of `v` must be promotable to those
of `A`. If so, yield the same result as `diagm(v) .+ A`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> addtodiag(A, v)
3×3 Array{Float64,2}:
 2.0  2.0   3.0
 2.0  6.0   6.0
 3.0  6.0  12.0
```

See also: [`addtodiag!`](@ref), [`subtractfromdiag!`](@ref), [`subtractfromdiag`](@ref)

"""
function addtodiag(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1})
    B = deepcopy(A);
    addtodiag!(B, v);
    return B;
end

"""

    subtractfromdiag!(A, v)

Efficiently subtract the main diagonal of the matrix `A` from the vector `v` in-place
(overwriting the latter). `length(v)` must match one of the dimensions of `A`, with the
other dimension being at least as large, and the element types of `v` must be promotable to
those of `A`. If so, yield the same result as `A .= diagm(v) .- A`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> subtractfromdiag!(A, v)
3×3 Array{Float64,2}:
  0.0  -2.0  -3.0
 -2.0  -2.0  -6.0
 -3.0  -6.0  -6.0
```

See also: [`addtodiag!`](@ref), [`addtodiag`](@ref), [`subtractfromdiag`](@ref)

"""
function subtractfromdiag!(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1})
    lmul!(-1.0, A);
    addtodiag!(A, v);
end


"""

    subtractfromdiag(A, v)

Efficiently subtract the main diagonal of the matrix `A` from the vector `v` not in-place,
returning the result. `length(v)` must match one of the dimensions of `A`, with the other
dimension being at least as large, and the element types of `v` must be promotable to those
of `A`. If so, yield the same result as `diagm(v) .- A`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> subtractfromdiag(A, v)
3×3 Array{Float64,2}:
 1.0  2.0  3.0
 2.0  4.0  6.0
 3.0  6.0  9.0
```

See also: [`addtodiag!`](@ref), [`addtodiag`](@ref), [`subtractfromdiag!`](@ref)

"""
function subtractfromdiag(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1})
    B = deepcopy(A);
    subtractfromdiag!(B, v);
    return B;
end

"""

    mulMatDiagleft!(A, v)

Efficiently multiply each row of the matrix `A` by the corresponding element of the vector
`v` in-place (overwriting the former). The element types of `v` must be promotable to those
of `A`. If so, yield the same result as `A .= diagm(v) * A`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> mulMatDiagleft!(A, v)
3×3 Array{Float64,2}:
 1.0   2.0   3.0
 4.0   8.0  12.0
 9.0  18.0  27.0
```

See also: [`mulMatDiagleft`](@ref), [`mulMatDiagright!`](@ref), [`mulMatDiagright`](@ref)

"""
mulMatDiagleft!(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1}) = broadcast!(*, A, A, v);

"""

    mulMatDiagleft(A, v)

Efficiently multiply each row of the matrix `A` by the corresponding element of the vector
`v` not in-place, returning the result. The element types of `v` must be promotable to those
of `A`. If so, yield the same result as `diagm(v) * A`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> mulMatDiagleft(A, v)
3×3 Array{Float64,2}:
 1.0   2.0   3.0
 4.0   8.0  12.0
 9.0  18.0  27.0
```

See also: [`mulMatDiagleft!`](@ref), [`mulMatDiagright!`](@ref), [`mulMatDiagright`](@ref)

"""
mulMatDiagleft(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1}) = broadcast(*, A, v);

"""

    mulMatDiagright!(A, v)

Efficiently multiply each column of the matrix `A` by the corresponding element of the
vector `v` in-place (overwriting the former). The element types of `v` must be promotable to
those of `A`, and `v` must have the same length as the number of columns of `A`. If so,
yield the same result as `A .= A * diagm(v)`.

# Examples
```jldoctest

julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> mulMatDiagright!(A, v)
3×3 Array{Float64,2}:
 1.0   4.0   9.0
 2.0   8.0  18.0
 3.0  12.0  27.0
```

See also: [`mulMatDiagleft!`](@ref), [`mulMatDiagleft`](@ref), [`mulMatDiagright`](@ref)

"""
function mulMatDiagright!(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1})
  
    if (length(v) != size(A, 2))
        error("mulMatDiagright!(): v (currently of length ", length(v),
              ") must have the same length as the number of columns in A (currently ",
              size(A, 2), ")");
    end
    
    @inbounds for jj=1:size(A, 2)
        @inbounds for ii=1:size(A, 1)
            A[ii, jj] *= v[jj];
        end
    end
    
end

"""

    mulMatDiagright(A, v)

Efficiently multiply each column of the matrix `A` by the corresponding element of the
vector `v` not in-place, returning the result. The element types of `v` must be promotable
to those of `A`. If so, yield the same result as `A * diagm(v)`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 2.0 4.0 6.0; 3.0 6.0 9.0];

julia> v = [1.0; 2.0; 3.0];

julia> mulMatDiagright(A, v)
3×3 Array{Float64,2}:
 1.0   4.0   9.0
 2.0   8.0  18.0
 3.0  12.0  27.0
```

See also: [`mulMatDiagleft!`](@ref), [`mulMatDiagleft`](@ref), [`mulMatDiagright!`](@ref)

"""
function mulMatDiagright(A::AbstractArray{<:Number, 2}, v::AbstractArray{<:Number, 1})
    B = deepcopy(A);
    mulMatDiagright!(B, v);
    return B;
end

"""

    vecNto3N(v)

Return a vector `w` that is 3 times the length of `v`, such that each 3-element block
`w[3*ii-2:3*ii]` has all elements identical to `v[ii]`. Mathematically, the vector
```math
v = \\begin{bmatrix}
v_1 &
v_2 &
\\ldots &
v_N
\\end{bmatrix}^{\\top}
```
is mapped to the vector
```math
w = \\begin{bmatrix}
v_1 &
v_1 &
v_1 &
v_2 &
v_2 &
v_2 &
\\ldots &
v_N &
v_N &
v_N
\\end{bmatrix}^{\\top}
```
# Examples
```jldoctest
julia> v = [1.0; 2.0; 3.0];

julia> vecNto3N(v)
9-element Array{Float64,1}:
 1.0
 1.0
 1.0
 2.0
 2.0
 2.0
 3.0
 3.0
 3.0
```

"""
vecNto3N(v::AbstractArray{<:Number, 1}) = vec(repeat(transpose(v), 3, 1));

"""

    sympart!(A)

Convert a square matrix `A` into its symmetric part in-place (overwriting `A`), yielding the
same result as `A .= (A .+ transpose(A))./2`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];

julia> B = 1im .* A;

julia> sympart!(A)
3×3 Array{Float64,2}:
 1.0  3.0  5.0
 3.0  5.0  7.0
 5.0  7.0  9.0

julia> sympart!(B)
3×3 Array{Complex{Float64},2}:
 0.0+1.0im  0.0+3.0im  0.0+5.0im
 0.0+3.0im  0.0+5.0im  0.0+7.0im
 0.0+5.0im  0.0+7.0im  0.0+9.0im
```

See also: [`sympart`](@ref), [`hermpart!`](@ref), [`hermpart`](@ref), [`ahermpart!`](@ref),
[`ahermpart`](@ref)

"""
function sympart!(A::AbstractArray{<:Number, 2})

    if (size(A, 1) != size(A, 2))
        error("sympart!(): A must be a square matrix");
    end
  
    @inbounds for jj=1:size(A, 2)
        @inbounds for ii=(jj + 1):size(A, 1)
            A[ii, jj] = (A[ii, jj] + A[jj, ii])/2;
            A[jj, ii] = A[ii, jj];
        end
    end

end

"""

    sympart(A)

Convert a square matrix `A` into its symmetric part not in-place, returning the result, and
yielding the same result as `(A .+ transpose(A))./2`.

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];

julia> B = 1im .* A;

julia> sympart(A)
3×3 Array{Float64,2}:
 1.0  3.0  5.0
 3.0  5.0  7.0
 5.0  7.0  9.0

julia> sympart(B)
3×3 Array{Complex{Float64},2}:
 0.0+1.0im  0.0+3.0im  0.0+5.0im
 0.0+3.0im  0.0+5.0im  0.0+7.0im
 0.0+5.0im  0.0+7.0im  0.0+9.0im
```

See also: [`sympart!`](@ref), [`hermpart!`](@ref), [`hermpart`](@ref), [`ahermpart!`](@ref),
[`ahermpart`](@ref)

"""
function sympart(A::AbstractArray{<:Number, 2})
    B = deepcopy(A);
    sympart!(B);
    return B;
end


"""

    hermpart!(A)

Convert a square matrix `A` into its Hermitian part in-place (overwriting `A`), yielding the
same result as `A .= (A .+ adjoint(A))./2`. (This automatically uses `sympart!()` if `A` is
real to avoid conflict with promotion to complex types.)

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];

julia> B = 1im .* A;

julia> hermpart!(A)
3×3 Array{Float64,2}:
 1.0  3.0  5.0
 3.0  5.0  7.0
 5.0  7.0  9.0

julia> hermpart!(B)
3×3 Array{Complex{Float64},2}:
 0.0+0.0im  0.0-1.0im  0.0-2.0im
 0.0+1.0im  0.0+0.0im  0.0-1.0im
 0.0+2.0im  0.0+1.0im  0.0+0.0im
```

See also: [`sympart!`](@ref), [`sympart`](@ref), [`hermpart`](@ref), [`ahermpart!`](@ref),
[`ahermpart`](@ref)

"""
function hermpart!(A::AbstractArray{<:Complex{<:Real}, 2})

    if (size(A, 1) != size(A, 2))
        error("hermpart!(): A must be a square matrix");
    end
  
    @inbounds for jj=1:size(A, 2)
        A[jj, jj] = real(A[jj, jj]);
        @inbounds for ii=(jj + 1):size(A, 1)
            A[ii, jj] = (A[ii, jj] + conj(A[jj, ii]))/2;
            A[jj, ii] = conj(A[ii, jj]);
        end
    end

end
hermpart!(A::AbstractArray{<:Real, 2}) = sympart!(A);

"""

    hermpart(A)

Convert a square matrix `A` into its Hermitian part not in-place, returning the result, and
yielding the same result as `(A .+ adjoint(A))./2`. (Automatically use `sympart()` if `A` is
real.)

# Examples
```jldoctest
julia> A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];

julia> B = 1im .* A;

julia> hermpart(A)
3×3 Array{Float64,2}:
 1.0  3.0  5.0
 3.0  5.0  7.0
 5.0  7.0  9.0

julia> hermpart(B)
3×3 Array{Complex{Float64},2}:
 0.0+0.0im  0.0-1.0im  0.0-2.0im
 0.0+1.0im  0.0+0.0im  0.0-1.0im
 0.0+2.0im  0.0+1.0im  0.0+0.0im
```

See also: [`sympart!`](@ref), [`sympart`](@ref), [`hermpart!`](@ref), [`ahermpart!`](@ref),
[`ahermpart`](@ref)

"""
function hermpart(A::AbstractArray{<:Complex{<:Real}, 2})
    B = deepcopy(A);
    hermpart!(B);
    return B;
end
hermpart(A::AbstractArray{<:Real, 2}) = sympart(A);

"""

    ahermpart!(A)

Convert a square matrix `A` into its antiHermitian part in-place (overwriting `A`), yielding
the same result as `A .= (A .- adjoint(A))./2im`. (Note that the result is Hermitian.)

# Examples
```jldoctest
julia> B = 1im .* [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];

julia> ahermpart!(B)
3×3 Array{Complex{Float64},2}:
 1.0+0.0im  3.0-0.0im  5.0-0.0im
 3.0+0.0im  5.0+0.0im  7.0-0.0im
 5.0+0.0im  7.0+0.0im  9.0+0.0im
```

See also: [`sympart!`](@ref), [`sympart`](@ref), [`hermpart!`](@ref), [`hermpart`](@ref),
[`ahermpart`](@ref)

"""
function ahermpart!(A::AbstractArray{<:Complex{<:Real}, 2})
    lmul!(-1im, A);
    hermpart!(A);
end

"""

    ahermpart(A)

Convert a square matrix `A` into its antiHermitian part not in-place, returning the result,
and yielding the same result as `(A .- adjoint(A))./2im`. (Note that the result is
Hermitian.)

# Examples
```jldoctest
julia> B = 1im .* [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];

julia> ahermpart(B)
3×3 Array{Complex{Float64},2}:
 1.0+0.0im  3.0-0.0im  5.0-0.0im
 3.0+0.0im  5.0+0.0im  7.0-0.0im
 5.0+0.0im  7.0+0.0im  9.0+0.0im
```

See also: [`sympart!`](@ref), [`sympart`](@ref), [`hermpart!`](@ref), [`hermpart`](@ref),
[`ahermpart!`](@ref)

"""
function ahermpart(A::AbstractArray{<:Complex{<:Real}, 2})
    B = deepcopy(A);
    ahermpart!(B);
    return B;
end

