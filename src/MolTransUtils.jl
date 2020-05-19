"""

    findCOM(pos0, startatom=1, endatom=size(pos0, 2))

Return the unweighted center of mass of a block `pos0[:, startatom:endatom]` of coordinates,
where each column `pos0[:, p]` has 3 elements representing the Cartesian coordinates of a
particle labeled ``p``.

"""
function findCOM(pos0::AbstractArray{<:AbstractFloat, 2}, startatom::Integer=1,
                 endatom::Integer=size(pos0, 2))
  
    if (size(pos0, 1) != 3)
        error("findCOM(): pos0 must by a 3-by-N real matrix");
    end
    if (startatom < 1)
        error("findCOM(): startatom must be >= 1");
    end
    if (endatom > size(pos0, 2))
        error("findCOM(): endatom must be <= N = size(pos0, 2)");
    end
    if (startatom > endatom)
        error("findCOM(): startatom must be <= endatom");
    end
    currCOM = zeros(eltype(pos0), 3);
    @inbounds for pp=startatom:endatom
        @inbounds for jj=1:3
            currCOM[jj] += pos0[jj, pp];
        end
    end

    return currCOM./(endatom - startatom + 1);

end

"""

    translateCOM!(pos0, dispCOM, startatom=1, endatom=size(pos0, 2))

Translate the block `pos0[:, startatom:endatom]` of coordinates in-place (overwriting it),
where each column `pos0[:, p]` has 3 elements representing the Cartesian coordinates of a
particle labeled ``p`` and `dispCOM` is a 3-element vector representing the displacement,
such that for each `p` in the range, `pos0[:, p] .+= dispCOM`.

See also: [`translateCOM`](@ref)

"""
function translateCOM!(pos0::AbstractArray{FT, 2}, dispCOM::Array{FT, 1},
                       startatom::Integer=1,
                       endatom::Integer=size(pos0, 2)) where FT<:AbstractFloat
  
    if (length(dispCOM) != 3)
        error("translateCOM!(): dispCOM must be a 3-element real vector");
    end
    if (size(pos0, 1) != 3)
        error("translateCOM!(): pos0 must by a 3-by-N real matrix");
    end
    if (startatom < 1)
        error("translateCOM!(): startatom must be >= 1");
    end
    if (endatom > size(pos0, 2))
        error("translateCOM!(): endatom must be <= N = size(pos0, 2)");
    end
    if (startatom > endatom)
        error("translateCOM!(): startatom must be <= endatom");
    end
  
    if (dispCOM != zero(dispCOM))
        @inbounds for pp=startatom:endatom
            @inbounds for jj=1:3
                pos0[jj, pp] += dispCOM[jj];
            end
        end
    end

end

"""

    rotatecoords!(pos0, rotmat, startatom=1, endatom=size(pos0, 2))

Rotate the block `pos0[:, startatom:endatom]` of coordinates in-place (overwriting it),
where each column `pos0[:, p]` has 3 elements representing the Cartesian coordinates of a
particle labeled ``p`` and `rotmat` is a 3-by-3 matrix representing the (alibi) rotation,
such that for each `p` in the range, `pos0[:, p] .= rotmat * pos0[:, p]`.

See also: [`rotatecoords`](@ref)

!!! warning
    
    This function will not check whether `rotmat` is real and orthogonal. It is the
    responsibility of the caller to ensure that `rotmat` represents a physical rotation;
    otherwise, unexpected behavior could emerge.

"""
function rotatecoords!(pos0::AbstractArray{FT, 2}, rotmat::AbstractArray{FT, 2},
                       startatom::Integer=1,
                       endatom::Integer=size(pos0, 2)) where FT<:AbstractFloat
    if (size(rotmat, 1) != 3 || size(rotmat, 2) != 3)
        error("rotatecoords!(): rotmat must be a 3-by-3 orthogonal real matrix");
    end
    if (size(pos0, 1) != 3)
        error("rotatecoords!(): pos0 must by a 3-by-N real matrix");
    end
    if (startatom < 1)
        error("rotatecoords!(): startatom must be >= 1");
    end
    if (endatom > size(pos0, 2))
        error("rotatecoords!(): endatom must be <= N = size(pos0, 2)");
    end
    if (startatom > endatom)
        error("rotatecoords!(): startatom must be <= endatom");
    end

    if (rotmat != Matrix{FT}(I, 3, 3))
        newvec = zero(pos0[:, 1]);
        @inbounds for pp=startatom:endatom
            setindex!(newvec, rotmat * pos0[:, pp], 1:3);
            @inbounds for jj=1:3
                pos0[jj, pp] = newvec[jj];
            end
        end
    end
    
end

"""

    translateCOM(pos0, dispCOM, startatom=1, endatom=size(pos0, 2))

Translate the block `pos0[:, startatom:endatom]` of coordinates not in-place, returning the
result, where each column `pos0[:, p]` has 3 elements representing the Cartesian coordinates
of a particle labeled ``p`` and `dispCOM` is a 3-element vector representing the
displacement, such that for each `p` in the range, `pos0[:, p] .+= dispCOM`.

See also: [`translateCOM!`](@ref)

"""
function translateCOM(pos0::AbstractArray{FT, 2}, dispCOM::Array{FT, 1},
                      startatom::Integer=1,
                      endatom::Integer=size(pos0, 2)) where FT<:AbstractFloat
    B = deepcopy(pos0);
    translateCOM!(B, dispCOM, startatom, endatom);
    return B;
end

"""

    rotatecoords(pos0, rotmat, startatom=1, endatom=size(pos0, 2))

Rotate the block `pos0[:, startatom:endatom]` of coordinates not in-place, returning the
result, where each column `pos0[:, p]` has 3 elements representing the Cartesian coordinates
of a particle labeled ``p`` and `rotmat` is a 3-by-3 matrix representing the rotation,
such that for each `p` in the range, `pos0[:, p] .= rotmat * pos0[:, p]`.

!!! warning

    This function will not check whether `rotmat` is real and orthogonal. It is the
    responsibility of the caller to ensure that `rotmat` represents a physical rotation;
    otherwise, unexpected behavior could emerge.

See also: [`rotatecoords!`](@ref)

"""
function rotatecoords(pos0::AbstractArray{FT, 2}, rotmat::AbstractArray{FT, 2},
                      startatom::Integer=1,
                      endatom::Integer=size(pos0, 2)) where FT<:AbstractFloat
    B = deepcopy(pos0);
    rotatecoords!(B, rotmat, startatom, endatom);
    return B;
end

"""

    rotate3Nmat!(A, rotmat, startatom=1, endatom=div(size(A, 2), 3))

Rotate the contents in the diagonal block of the complex matrix
`A[3*startatom-2:3*endatom, 3*startatom-2:3*endatom]` by rotating each 3-by-3 block therein
by the 3-by-3 real orthogonal (alibi) rotation matrix `rotmat`; this is done in-place,
overwriting `A`.

See also: [`rotate3Nmat`](@ref)

!!! warning
    
    This function will not check whether `rotmat` is real and orthogonal. It is the
    responsibility of the caller to ensure that `rotmat` represents a physical rotation;
    otherwise, unexpected behavior could emerge.

"""
function rotate3Nmat!(A::AbstractArray{Complex{FT}, 2}, rotmat::AbstractArray{FT, 2},
                      startatom::Integer=1,
                      endatom::Integer=div(size(A, 2), 3)) where FT<:AbstractFloat
  
    if (size(rotmat, 1) != 3 || size(rotmat, 2) != 3)
        error("rotate3Nmat!(): rotmat must be a 3-by-3 orthogonal real matrix");
    end
    if (size(A, 1) != size(A, 2) || mod(size(A, 1), 3) != 0)
        error("rotate3Nmat!(): A must by a 3N-by-3N real or complex matrix");
    end
    if (startatom < 1)
        error("rotate3Nmat!(): startatom must be >= 1");
    end
    if (endatom > div(size(A, 2), 3))
        error("rotate3Nmat!(): endatom must be <= N = size(A, 2)/3");
    end
    if (startatom > endatom)
        error("rotate3Nmat!(): startatom must be <= endatom");
    end

    if (rotmat != Matrix{FT}(I, 3, 3))
        @inbounds for qq=startatom:endatom
            @inbounds for pp=startatom:endatom
                tempmat = view(A, 3*pp-2:3*pp, 3*qq-2:3*qq);
                tempmat .= rotmat * tempmat * transpose(rotmat);
            end
        end
    end

end

"""

    rotate3Nmat(A, rotmat, startatom=1, endatom=div(size(A, 2), 3))

Rotate the contents in the diagonal block of the complex matrix
`A[3*startatom-2:3*endatom, 3*startatom-2:3*endatom]` by rotating each 3-by-3 block therein
by the 3-by-3 real orthogonal (alibi) rotation matrix `rotmat`; this is done not in-place,
returning the result.

See also: [`rotate3Nmat!`](@ref)

!!! warning
    
    This function will not check whether `rotmat` is real and orthogonal. It is the
    responsibility of the caller to ensure that `rotmat` represents a physical rotation;
    otherwise, unexpected behavior could emerge.

"""
function rotate3Nmat(A::AbstractArray{Complex{FT}, 2}, rotmat::AbstractArray{FT, 2},
                     startatom::Integer=1,
                     endatom::Integer=div(size(A, 2), 3)) where FT<:AbstractFloat
    B = deepcopy(A);
    rotate3Nmat!(B, rotmat, startatom, endatom);
    return B;
end

"""

    makeCOM0!(pos0, startatom=1, endatom=size(pos0, 2))

Translate the block `pos0[:, startatom:endatom]` of coordinates in-place (overwriting it) to
have a center of mass of 0, where each column `pos0[:, p]` has 3 elements representing the
Cartesian coordinates of a particle labeled ``p``.

See also: [`makeCOM0`](@ref)

"""
function makeCOM0!(pos0::AbstractArray{<:AbstractFloat, 2}, startatom::Integer=1,
                   endatom::Integer=size(pos0, 2))
  
    currCOM = findCOM(pos0, startatom, endatom);
    lmul!(-1, currCOM);
    translateCOM!(pos0, currCOM, startatom, endatom);
    
end

"""

    makeCOM0(pos0, startatom=1, endatom=size(pos0, 2))

Translate the block `pos0[:, startatom:endatom]` of coordinates not in-place, returning the
result, to have a center of mass of 0, where each column `pos0[:, p]` has 3 elements
representing the Cartesian coordinates of a particle labeled ``p``.

See also: [`makeCOM0!`](@ref)

"""
function makeCOM0(pos0::AbstractArray{<:AbstractFloat, 2}, startatom::Integer=1,
                  endatom::Integer=size(pos0, 2))

    currCOM = findCOM(pos0, startatom, endatom);
    lmul!(-1.0, currCOM);
    return translateCOM(pos0, currCOM, startatom, endatom);

end
