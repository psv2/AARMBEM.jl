## The first set of functions in this file are baseline functions that
##  modify a matrix representing the Green's function, operating at
##  complex frequency

"""

    GFPECGG!(GscaGG::AbstractArray{Complex{FT}, 2},
             startidx1::Integer, startidx2::Integer,
             freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
             posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
             latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
             alphafield::FT, alphasource::FT, numdims::Integer) where FT<:AbstractFloat

Fill `GscaGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 scattering
Green's function interaction tensor above a perfect electrically conducting (PEC) plane
(coplanar with the xy-plane) from the vacuum Green's function via image theory, at frequency
`freq` and wavevector `blochk` (a 3-element vector), corresponding to source position
`possource` and field position `posfield` (both 3-element vectors) described by Gaussian
basis functions parameterized by real atomic polarizabilities `alphasource` and `alphafield`
respectively, using the appropriate function for `numdims` periodic dimensions (with lattice
parameters `latticevecs` for the real lattice vectors and `reciprocalvecs` for the
reciprocal lattice vectors as appropriate). `numdims` must be 0, 1, or 2.

"""
function GFPECGG!(GscaGG::AbstractArray{Complex{FT}, 2},
                  startidx1::Integer, startidx2::Integer,
                  freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
                  alphafield::FT, alphasource::FT, numdims::Integer) where FT<:AbstractFloat
  
    if (numdims == 0)
        GFPECGG0!(GscaGG, startidx1, startidx2, freq, posfield, possource, alphafield, alphasource);
    elseif (numdims == 1)
        GFPECGG1Ewald!(GscaGG, startidx1, startidx2, freq, blochk, posfield, possource, latticevecs, reciprocalvecs, alphafield, alphasource);
    elseif (numdims == 2)
        GFPECGG2Ewald!(GscaGG, startidx1, startidx2, freq, blochk, posfield, possource, latticevecs, reciprocalvecs, alphafield, alphasource);
    else
        error("GFPECGG!(): Periodic dimensionality numdims (currently ", numdims,
              ") must be 0 (no periodicity), 1, or 2");
    end
    
end

"""

    GFPECGG0!(GscaGG::AbstractArray{Complex{FT}, 2},
              startidx1::Integer, startidx2::Integer,
              freq::Union{FT, Complex{FT}},
              posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
              alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Fill `GscaGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 scattering
Green's function interaction tensor above a perfect electrically conducting (PEC) plane
(coplanar with the xy-plane) from the vacuum Green's function via image theory, at frequency
`freq` without spatial periodicity, corresponding to source position `possource` and field
position `posfield` (both 3-element vectors) described by Gaussian basis functions
parameterized by real atomic polarizabilities `alphasource` and `alphafield` respectively.

"""
function GFPECGG0!(GscaGG::AbstractArray{Complex{FT}, 2},
                   startidx1::Integer, startidx2::Integer,
                   freq::Union{FT, Complex{FT}},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    if (posfield[3] < zero(eltype(posfield)))
        error("GFPECGG0!(): can only handle field positions z (currently ",
              posfield[3], ") larger than 0");
    end
    if (possource[3] < zero(eltype(possource)))
        error("GFPECGG0!(): can only handle source positions z' (currently ",
              possource[3], ") larger than 0");
    end

    possource[3] *= -1;
    GFVACGG0!(GscaGG, startidx1, startidx2, freq,
              posfield, possource, alphafield, alphasource);

    for jj=1:2
        for ii=1:3
            GscaGG[ii, jj] *= -1;
        end
    end

    possource[3] *= -1;

end

"""

    GFPECGG1Ewald!(GscaGG::AbstractArray{Complex{FT}, 2},
                   startidx1::Integer, startidx2::Integer,
                   freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   latticevecs::AbstractArray{FT, 2},
                   reciprocalvecs::AbstractArray{FT, 2},
                   alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Fill `GscaGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 scattering
Green's function interaction tensor above a perfect electrically conducting (PEC) plane
(coplanar with the xy-plane) from the vacuum Green's function via image theory, at frequency
`freq` and wavevector `blochk` (a 3-element vector), corresponding to source position
`possource` and field position `posfield` (both 3-element vectors) described by Gaussian
basis functions parameterized by real atomic polarizabilities `alphasource` and `alphafield`
respectively, in 1 periodic dimension (with lattice parameters `latticevecs` for the real
lattice vectors and `reciprocalvecs` for the reciprocal lattice vectors as appropriate).
This is done via Ewald summation over Gaussian basis functions in 1 dimension.

!!! warning
    
    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFPECGG1Ewald!(GscaGG::AbstractArray{Complex{FT}, 2},
                        startidx1::Integer, startidx2::Integer,
                        freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                        posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                        latticevecs::AbstractArray{FT, 2},
                        reciprocalvecs::AbstractArray{FT, 2},
                        alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    if (posfield[3] < zero(eltype(posfield)))
        error("GFPECGG1Ewald!(): can only handle field positions z (currently ",
              posfield[3], ") larger than 0");
    end
    if (possource[3] < zero(eltype(possource)))
        error("GFPECGG1Ewald!(): can only handle source positions z' (currently ",
              possource[3], ") larger than 0");
    end
    if (abs(latticevecs[3, 1])/norm(latticevecs[:, 1]) > eps())
        error("GFPECGG1Ewald!(): 1D-periodic lattice vector must be",
              " perpendicular to e_z = [0, 0, 1]");
    end

    possource[3] *= -1;
    GFVACGG1Ewald!(GscaGG, startidx1, startidx2, freq, blochk, posfield, possource,
                   latticevecs, reciprocalvecs, alphafield, alphasource);
    
    for jj=1:2
        for ii=1:3
            GscaGG[ii, jj] *= -1;
        end
    end
    
    possource[3] *= -1;

end

"""

    GFPECGG2Ewald!(GscaGG::AbstractArray{Complex{FT}, 2},
                   startidx1::Integer, startidx2::Integer,
                   freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   latticevecs::AbstractArray{FT, 2},
                   reciprocalvecs::AbstractArray{FT, 2},
                   alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Fill `GscaGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 scattering
Green's function interaction tensor above a perfect electrically conducting (PEC) plane
(coplanar with the xy-plane) from the vacuum Green's function via image theory, at frequency
`freq` and wavevector `blochk` (a 3-element vector), corresponding to source position
`possource` and field position `posfield` (both 3-element vectors) described by Gaussian
basis functions parameterized by real atomic polarizabilities `alphasource` and `alphafield`
respectively, in 2 periodic dimensions (with lattice parameters `latticevecs` for the real
lattice vectors and `reciprocalvecs` for the reciprocal lattice vectors as appropriate).
This is done via Ewald summation over Gaussian basis functions in 2 dimensions.

!!! warning
    
    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFPECGG2Ewald!(GscaGG::AbstractArray{Complex{FT}, 2},
                        startidx1::Integer, startidx2::Integer,
                        freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                        posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                        latticevecs::AbstractArray{FT, 2},
                        reciprocalvecs::AbstractArray{FT, 2},
                        alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    if (posfield[3] < zero(eltype(posfield)))
        error("GFPECGG2Ewald!(): can only handle field positions z (currently ",
              posfield[3], ") larger than 0");
    end
    if (possource[3] < zero(eltype(possource)))
        error("GFPECGG2Ewald!(): can only handle source positions z' (currently ",
              possource[3], ") larger than 0");
    end
    if (abs(latticevecs[3, 1])/norm(latticevecs[:, 1]) > eps() ||
        abs(latticevecs[3, 2])/norm(latticevecs[:, 2]) > eps())
        error("GFPECGG2Ewald!(): 2D-periodic lattice vectors must be perpendicular",
              " to e_z = [0, 0, 1]");
    end

    possource[3] *= -1;
    GFVACGG2Ewald!(GscaGG, startidx1, startidx2, freq, blochk, posfield, possource,
                   latticevecs, reciprocalvecs, alphafield, alphasource);

    for jj=1:2
        for ii=1:3
            GscaGG[ii, jj] *= -1;
        end
    end
    
    possource[3] *= -1;

end

## The second set of functions in this file consists of aliases for
##  functions with fewer arguments, useful for compact geometries

"""

    GFPECGG!(GscaGG::AbstractArray{Complex{FT}, 2},
             startidx1::Integer, startidx2::Integer,
             freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
             possource::AbstractArray{FT, 1}, alphafield::FT,
             alphasource::FT) where {FT<:AbstractFloat}

Fill `GscaGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 scattering
Green's function interaction tensor above a perfect electrically conducting (PEC) plane
(coplanar with the xy-plane) from the vacuum Green's function via image theory, at frequency
`freq`, corresponding to source position `possource` and field position `posfield` (both
3-element vectors) described by Gaussian basis functions parameterized by real atomic
polarizabilities `alphasource` and `alphafield` respectively, for nonperiodic geometries.

"""
GFPECGG!(GscaGG::AbstractArray{Complex{FT}, 2},
         startidx1::Integer, startidx2::Integer,
         freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
         possource::AbstractArray{FT, 1}, alphafield::FT,
         alphasource::FT) where {FT<:AbstractFloat} =
GFPECGG0!(GscaGG, startidx1, startidx2, freq, posfield, possource, alphafield, alphasource);

## The third set of functions in this file consists of aliases for not
##  in-place function versions of the first and second sets
## (Note that all return types are 3-by-3 arrays with Complex{Float64}
##  elements, independent of arguments)

"""

    GFPECGG(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                 posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                 latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
                 alphafield::FT, alphasource::FT, numdims::Integer) where FT<:AbstractFloat

Return the 3-by-3 scattering Green's function interaction tensor above a perfect
electrically conducting (PEC) plane (coplanar with the xy-plane) from the vacuum Green's
function via image theory, at frequency `freq` and wavevector `blochk` (a 3-element vector),
corresponding to source position `possource` and field position `posfield` (both 3-element
vectors) described by Gaussian basis functions parameterized by real atomic polarizabilities
`alphasource` and `alphafield` respectively, using the appropriate function for `numdims`
periodic dimensions (with lattice parameters `latticevecs` for the real lattice vectors and
`reciprocalvecs` for the reciprocal lattice vectors as appropriate). `numdims` must be 0, 1,
or 2.

"""
function GFPECGG(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                 posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                 latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
                 alphafield::FT, alphasource::FT, numdims::Integer) where FT<:AbstractFloat
  
    GscaGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFPECGG!(GscaGG, 1, 1, freq, blochk, posfield, possource,
             latticevecs, reciprocalvecs, alphafield, alphasource, numdims);
    return GscaGG;
    
end

"""

    GFPECGG0(freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
             possource::AbstractArray{FT, 1}, alphafield::FT,
             alphasource::FT) where FT<:AbstractFloat

Return the 3-by-3 scattering Green's function interaction tensor above a perfect
electrically conducting (PEC) plane (coplanar with the xy-plane) from the vacuum Green's
function via image theory, at frequency `freq` without spatial periodicity, corresponding to
source position `possource` and field position `posfield` (both 3-element vectors) described
by Gaussian basis functions parameterized by real atomic polarizabilities `alphasource` and
`alphafield` respectively.

"""
function GFPECGG0(freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
                  possource::AbstractArray{FT, 1}, alphafield::FT,
                  alphasource::FT) where FT<:AbstractFloat
  
    GscaGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFPECGG0!(GscaGG, 1, 1, freq, posfield, possource, alphafield, alphasource);
    return GscaGG;
    
end

"""

    GFPECGG1Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  latticevecs::AbstractArray{FT, 2},
                  reciprocalvecs::AbstractArray{FT, 2}, alphafield::FT,
                  alphasource::FT) where FT<:AbstractFloat

Return the 3-by-3 scattering Green's function interaction tensor above a perfect
electrically conducting (PEC) plane (coplanar with the xy-plane) from the vacuum Green's
function via image theory, at frequency `freq` and wavevector `blochk` (a 3-element vector),
corresponding to source position `possource` and field position `posfield` (both 3-element
vectors) described by Gaussian basis functions parameterized by real atomic polarizabilities
`alphasource` and `alphafield` respectively, in 1 periodic dimension (with lattice parameters
`latticevecs` for the real lattice vectors and `reciprocalvecs` for the reciprocal lattice
vectors as appropriate). This is done via Ewald summation over Gaussian basis functions in 1
dimension.

!!! warning

    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFPECGG1Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                       posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                       latticevecs::AbstractArray{FT, 2},
                       reciprocalvecs::AbstractArray{FT, 2}, alphafield::FT,
                       alphasource::FT) where FT<:AbstractFloat

    GscaGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFPECGG1Ewald!(GscaGG, 1, 1, freq, blochk, posfield, possource,
                   latticevecs, reciprocalvecs, alphafield, alphasource);
    return GscaGG;

end

"""

    GFPECGG2Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  latticevecs::AbstractArray{FT, 2},
                  reciprocalvecs::AbstractArray{FT, 2}, alphafield::FT,
                  alphasource::FT) where FT<:AbstractFloat

Return the 3-by-3 scattering Green's function interaction tensor above a perfect
electrically conducting (PEC) plane (coplanar with the xy-plane) from the vacuum Green's
function via image theory, at frequency `freq` and wavevector `blochk` (a 3-element vector),
corresponding to source position `possource` and field position `posfield` (both 3-element
vectors) described by Gaussian basis functions parameterized by real atomic polarizabilities
`alphasource` and `alphafield` respectively, in 2 periodic dimensions (with lattice
parameters `latticevecs` for the real lattice vectors and `reciprocalvecs` for the reciprocal
lattice vectors as appropriate). This is done via Ewald summation over Gaussian basis
functions in 2 dimensions.

!!! warning
    
    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFPECGG2Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                       posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                       latticevecs::AbstractArray{FT, 2},
                       reciprocalvecs::AbstractArray{FT, 2}, alphafield::FT,
                       alphasource::FT) where FT<:AbstractFloat

    GscaGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFPECGG2Ewald!(GscaGG, 1, 1, freq, blochk, posfield, possource,
                   latticevecs, reciprocalvecs, alphafield, alphasource);
    return GscaGG;

end

"""

    GFPECGG(freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
            possource::AbstractArray{FT, 1}, alphafield::FT,
            alphasource::FT) where {FT<:AbstractFloat}

Return the 3-by-3 scattering Green's function interaction tensor above a perfect
electrically conducting (PEC) plane (coplanar with the xy-plane) from the vacuum Green's
function via image theory, at frequency `freq`, corresponding to source position `possource`
and field position `posfield` (both 3-element vectors) described by Gaussian basis functions
parameterized by real atomic polarizabilities `alphasource` and `alphafield` respectively,
for nonperiodic geometries.

"""
GFPECGG(freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
         possource::AbstractArray{FT, 1}, alphafield::FT,
         alphasource::FT) where {FT<:AbstractFloat} =
GFPECGG0(freq, posfield, possource, alphafield, alphasource);
