"""
The vacuum speed of light is a global constant (299792458 meters per second).
"""
const c = 299792458;

## The first set of functions in this file are baseline functions that
##  modify (or are required by those that modify) a matrix representing
##  the Green's function, operating at complex frequency

"""

    GFVACGG!(GvacGG::AbstractArray{Complex{FT}, 2},
             startidx1::Integer, startidx2::Integer,
             freq::Union{FT, Complex{FT}}, blochk::Array{FT, 1},
             posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
             latticevecs::Array{FT, 2}, reciprocalvecs::Array{FT, 2},
             alphafield::FT, alphasource::FT,
             numdims::Integer) where FT<:AbstractFloat

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` and wavevector `blochk` (a
3-element vector), corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively, using the
appropriate function for `numdims` periodic dimensions (with lattice parameters
`latticevecs` for the real lattice vectors and `reciprocalvecs` for the reciprocal lattice
vectors as appropriate). `numdims` must be 0, 1, or 2.

"""
function GFVACGG!(GvacGG::AbstractArray{Complex{FT}, 2},
                  startidx1::Integer, startidx2::Integer,
                  freq::Union{FT, Complex{FT}}, blochk::Array{FT, 1},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  latticevecs::Array{FT, 2}, reciprocalvecs::Array{FT, 2},
                  alphafield::FT, alphasource::FT,
                  numdims::Integer) where FT<:AbstractFloat

    if (numdims == 0)
        GFVACGG0!(GvacGG, startidx1, startidx2, freq,
                  posfield, possource, alphafield, alphasource);
    elseif (numdims == 1)
        GFVACGG1Ewald!(GvacGG, startidx1, startidx2, freq, blochk,
                       posfield, possource, latticevecs, reciprocalvecs,
                       alphafield, alphasource);
    elseif (numdims == 2)
        GFVACGG2Ewald!(GvacGG, startidx1, startidx2, freq, blochk,
                       posfield, possource, latticevecs, reciprocalvecs,
                       alphafield, alphasource);
    else
        error("GFVACGG!(): Periodic dimensionality numdims (currently ", numdims,
              ") must be 0 (no periodicity), 1, or 2 for GFVACGG");
    end

end

"""

    GFVACGG_sca!(GvacGG::AbstractArray{Complex{FT}, 2},
                 startidx1::Integer, startidx2::Integer,
                 freq::Union{FT, Complex{FT}}, blochk::Array{FT, 1},
                 posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                 latticevecs::Array{FT, 2}, reciprocalvecs::Array{FT, 2},
                 alphafield::FT, alphasource::FT,
                 numdims::Integer) where FT <: AbstractFloat

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with a 3-by-3 block of
zeros, which is the scattering Green's function in vacuum by definition, irrespective of the
other function arguments. 

"""
function GFVACGG_sca!(GvacGG::AbstractArray{Complex{FT}, 2},
                      startidx1::Integer, startidx2::Integer,
                      freq::Union{FT, Complex{FT}}, blochk::Array{FT, 1},
                      posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                      latticevecs::Array{FT, 2}, reciprocalvecs::Array{FT, 2},
                      alphafield::FT, alphasource::FT,
                      numdims::Integer) where FT <: AbstractFloat
  
    if (length(size(GvacGG)) != 2 || size(GvacGG, 1) < 3 || size(GvacGG, 2) < 3)
        error("GFVACGG_sca!(): GvacGG (current size: ",
              size(GvacGG), ") must be 3-by-3");
    end
    if (startidx1 < 1 || startidx1 > size(GvacGG, 1) - 2)
        error("GFVACGG_sca!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(GvacGG, 1) - 2,
              " as GvacGG has ", size(GvacGG, 1), " rows");
    end
    if (startidx2 < 1 || startidx2 > size(GvacGG, 2) - 2)
        error("GFVACGG_sca!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(GvacGG, 2) - 2,
              " as GvacGG has ", size(GvacGG, 2), " columns");
    end
    GvacGGtemp = view(GvacGG, startidx1:startidx1+2, startidx2:startidx2+2);
    fill!(GvacGGtemp, zero(eltype(GvacGGtemp)));

end

"""

    GFVACGG0!(GvacGG::AbstractArray{Complex{FT}, 2},
              startidx1::Integer, startidx2::Integer,
              freq::Union{FT, Complex{FT}},
              posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
              alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Fill `GvacGG[startidx1:startidx1+2,startidx2: startidx2+2]` in-place with the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` without spatial periodicity,
corresponding to source position `possource` and field position `posfield` (both 3-element
vectors) described by Gaussian basis functions parameterized by real atomic polarizabilities
`alphasource` and `alphafield` respectively.

"""
function GFVACGG0!(GvacGG::AbstractArray{Complex{FT}, 2},
                   startidx1::Integer, startidx2::Integer,
                   freq::Union{FT, Complex{FT}},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    if (length(size(GvacGG)) != 2 || size(GvacGG, 1) < 3 || size(GvacGG, 2) < 3)
        error("GFVACGG0!(): GvacGG (current size: ", size(GvacGG), ") must be 3-by-3");
    end
    if (startidx1 < 1 || startidx1 > size(GvacGG, 1) - 2)
        error("GFVACGG0!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(GvacGG, 1) - 2,
              " as GvacGG has ", size(GvacGG, 1), " rows");
    end
    if (startidx2 < 1 || startidx2 > size(GvacGG, 2) - 2)
        error("GFVACGG0!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(GvacGG, 2) - 2,
              " as GvacGG has ", size(GvacGG, 2), " columns");
    end
  
  # sigma_p = (alpha_p / 3)^(1/3) / (2*sqrt(pi)), sigmatot = sqrt(2(sigma_p^2 + sigma_q^2))
    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));

    myrhat = posfield .- possource;
    rho = norm(myrhat)/sigmatot;
    GvacGGtemp = view(GvacGG, startidx1:startidx1+2, startidx2:startidx2+2);
    fill!(GvacGGtemp, zero(eltype(GvacGGtemp)));
    if (rho > eps())
      
        q = sigmatot*freq/c;
        normalize!(myrhat);
        
        term1 = exp(-1*(q^2 / 4) + 1im*rho*q) * erfc(-1im*(q/2) - rho);
        if (isnan(term1)) # just in case
            term1 = zero(eltype(GvacGGtemp));
        end
# exp(-0.25*q^2 - 1im*q*rho)*erfc(-0.5im*q + rho), computed with better numerical stability given rho >= 0
        term2 = exp(-1*(q^2 / 4) - 1im*q*rho - (-1im*(q/2) +
                                                rho)^2 + log(erfcx(-1im*(q/2) + rho)));
        if (isnan(term2)) # just in case
            term2 = zero(eltype(GvacGGtemp));
        end

        gvacGG = (term1 - term2)/(8*pi*sigmatot*rho);
        dgvacGGdrho = ((4/sqrt(pi))*rho*exp(-1*rho^2) +
                       (-1 + 1im*q*rho)*term1 +
                       (1 + 1im*q*rho)*term2)/(8*pi*sigmatot*rho^2);
        d2gvacGGdrho2 = ((-8/sqrt(pi))*(rho + rho^3)*exp(-1*rho^2) -
                         (-2 + 2im*q*rho + (q*rho)^2)*term1 +
                         (-2 - 2im*q*rho + (q*rho)^2)*term2)/(8*pi*sigmatot*rho^3);

        GvacGGtemp .+= (d2gvacGGdrho2/(sigmatot^2) -
                        dgvacGGdrho/(sigmatot^2 * rho)).*complex(myrhat*transpose(myrhat)) .+
        Matrix{Complex{FT}}((dgvacGGdrho/(sigmatot^2 * rho) + (freq/c)^2 * gvacGG)*I, 3, 3);
    
        if (sum(isnan.(GvacGGtemp)) > 0)
            fill!(GvacGGtemp, zero(eltype(GvacGGtemp)));
        end
        
    end
  
end

"""

    GFVACGG0_sca!(GvacGG::AbstractArray{Complex{FT}, 2},
                  startidx1::Integer, startidx2::Integer,
                  freq::Union{FT, Complex{FT}},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  alphafield::FT, alphasource::FT) where FT <: AbstractFloat

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with a 3-by-3 block of
zeros, which is the scattering Green's function in vacuum by definition, irrespective of the
other function arguments. 

"""
function GFVACGG0_sca!(GvacGG::AbstractArray{Complex{FT}, 2},
                       startidx1::Integer, startidx2::Integer,
                       freq::Union{FT, Complex{FT}},
                       posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                       alphafield::FT, alphasource::FT) where FT<:AbstractFloat
  
    if (length(size(GvacGG)) != 2 || size(GvacGG, 1) < 3 || size(GvacGG, 2) < 3)
        error("GFVACGG0_sca!(): GvacGG (current size: ", size(GvacGG), ") must be 3-by-3");
    end
    if (startidx1 < 1 || startidx1 > size(GvacGG, 1) - 2)
        error("GFVACGG0_sca!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(GvacGG, 1) - 2,
              " as GvacGG has ", size(GvacGG, 1), " rows");
    end
    if (startidx2 < 1 || startidx2 > size(GvacGG, 2) - 2)
        error("GFVACGG: startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(GvacGG, 2) - 2,
              " as GvacGG has ", size(GvacGG, 2), " columns");
    end
    GvacGGtemp = view(GvacGG, startidx1:startidx1+2, startidx2:startidx2+2);
    fill!(GvacGGtemp, zero(eltype(GvacGGtemp)));
    
end

"""

    GFVACGGcoincident!(GvacGG::AbstractArray{Complex{FT}, 2},
                       startidx1::Integer, startidx2::Integer,
                       freq::Union{FT, Complex{FT}}, alphafield::FT,
                       alphasource::FT) where FT<:AbstractFloat

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 tensor
block corresponding to vacuum Green's function interaction between two Gaussian basis
functions, defined by respective atomic polarizabilities `alphasource` and `alphafield`,
when their centers coincide; this quantity always has a finite imaginary part, and has a
finite real part when `alphafield` and `alphasource` are not both zero (meaning at least one
basis function has a finite width).

"""
function GFVACGGcoincident!(GvacGG::AbstractArray{Complex{FT}, 2},
                            startidx1::Integer, startidx2::Integer,
                            freq::Union{FT, Complex{FT}}, alphafield::FT,
                            alphasource::FT) where FT<:AbstractFloat
  
    if (length(size(GvacGG)) != 2 || size(GvacGG, 1) < 3 || size(GvacGG, 2) < 3)
        error("GFVACGGcoincident!(): GvacGG (current size: ",
              size(GvacGG), ") must be 3-by-3");
    end
    if (startidx1 < 1 || startidx1 > size(GvacGG, 1) - 2)
        error("GFVACGGcoincident!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(GvacGG, 1) - 2,
              " as GvacGG has ", size(GvacGG, 1), " rows");
    end
    if (startidx2 < 1 || startidx2 > size(GvacGG, 2) - 2)
        error("GFVACGGcoincident!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(GvacGG, 2) - 2,
              " as GvacGG has ", size(GvacGG, 2), " columns");
    end

    setindex!(GvacGG,
              Matrix{eltype(GvacGG)}(GFVACGGcoincident(freq, alphafield, alphasource)*I,
                                     3, 3), startidx1:startidx1+2, startidx2:startidx2+2);

end

"""

    GFVACGG1Ewald!(GvacGG::AbstractArray{Complex{FT}, 2},
                   startidx1::Integer, startidx2::Integer,
                   freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   latticevecs::AbstractArray{FT, 2},
                   reciprocalvecs::AbstractArray{FT, 2},
                   alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` and wavevector `blochk` (a
3-element vector), corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively, in 1 periodic
dimension (with lattice parameters `latticevecs` for the real lattice vectors and
`reciprocalvecs` for the reciprocal lattice vectors as appropriate). This is done via Ewald
summation over Gaussian basis functions in 1 dimension.

!!! warning
    
    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFVACGG1Ewald!(GvacGG::AbstractArray{Complex{FT}, 2},
                        startidx1::Integer, startidx2::Integer,
                        freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                        posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                        latticevecs::AbstractArray{FT, 2},
                        reciprocalvecs::AbstractArray{FT, 2},
                        alphafield::FT, alphasource::FT) where FT<:AbstractFloat
  
    if (length(size(GvacGG)) != 2 || size(GvacGG, 1) < 3 || size(GvacGG, 2) < 3)
        error("GFVACGG1Ewald!(): GvacGG (current size: ", size(GvacGG), ") must be 3-by-3");
    end
    if (startidx1 < 1 || startidx1 > size(GvacGG, 1) - 2)
        error("GFVACGG1Ewald!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(GvacGG, 1) - 2,
              " as GvacGG has ", size(GvacGG, 1), " rows");
    end
    if (startidx2 < 1 || startidx2 > size(GvacGG, 2) - 2)
        error("GFVACGG1Ewald!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(GvacGG, 2) - 2,
              " as GvacGG has ", size(GvacGG, 2), " columns");
    end

    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));
  
    myrvec = posfield .- possource;
    mylatticeunitvec = normalize(latticevecs[:,1]);
    myrperp = norm(myrvec .- mylatticeunitvec .* dot(mylatticeunitvec, myrvec));
    kappa = min(sqrt(pi)/norm(latticevecs[:,1]), 1/myrperp);
    N1max = round(Int, max(10, 10*myrperp/sigmatot));
    N2max = 0;
    Mmax = 10;
    Smax = Mmax;
    GvacGGtemp = view(GvacGG, startidx1:startidx1+2, startidx2:startidx2+2);
    fill!(GvacGGtemp, zero(Complex{FT}));
    
    if (norm(myrvec)/sigmatot <= eps())
        GvacGGtemp .+= GFVACGGEwaldSR(freq, blochk, posfield, possource,
                                      latticevecs, alphafield, alphasource,
                                      kappa, N1max, N2max) .+
        GFVACGGEwaldLR1D(freq, blochk, posfield, possource, reciprocalvecs[:,1],
                         alphafield, alphasource, kappa, Mmax, Smax) -
        GFVACGGcoincident(freq, alphafield, alphasource)*I;
    else
        GvacGGtemp .+= GFVACGGEwaldSR(freq, blochk, posfield, possource,
                                      latticevecs, alphafield, alphasource,
                                      kappa, N1max, N2max) .+
        GFVACGGEwaldLR1D(freq, blochk, posfield, possource, reciprocalvecs[:,1],
                         alphafield, alphasource, kappa, Mmax, Smax);
    end

end

"""

    GFVACGG2Ewald!(GvacGG::AbstractArray{Complex{FT}, 2},
                   startidx1::Integer, startidx2::Integer,
                   freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   latticevecs::AbstractArray{FT, 2},
                   reciprocalvecs::AbstractArray{FT, 2},
                   alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` and wavevector `blochk` (a
3-element vector), corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively, in 2 periodic
dimensions (with lattice parameters `latticevecs` for the real lattice vectors and
`reciprocalvecs` for the reciprocal lattice vectors as appropriate). This is done via Ewald
summation over Gaussian basis functions in 2 dimensions.

!!! warning
    
    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFVACGG2Ewald!(GvacGG::AbstractArray{Complex{FT}, 2},
                        startidx1::Integer, startidx2::Integer,
                        freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                        posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                        latticevecs::AbstractArray{FT, 2},
                        reciprocalvecs::AbstractArray{FT, 2},
                        alphafield::FT, alphasource::FT) where FT<:AbstractFloat
  
    if (length(size(GvacGG)) != 2 || size(GvacGG, 1) < 3 || size(GvacGG, 2) < 3)
        error("GFVACGG2Ewald!(): GvacGG (current size: ", size(GvacGG), ") must be 3-by-3");
    end
    if (startidx1 < 1 || startidx1 > size(GvacGG, 1) - 2)
        error("GFVACGG2Ewald!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ", size(GvacGG, 1) - 2,
              " as GvacGG has ", size(GvacGG, 1), " rows");
    end
    if (startidx2 < 1 || startidx2 > size(GvacGG, 2) - 2)
        error("GFVACGG2Ewald!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ", size(GvacGG, 2) - 2,
              " as GvacGG has ", size(GvacGG, 2), " columns");
    end

    Auc = norm(cross(latticevecs[:,1], latticevecs[:,2]));
    M1max = 5;
    M2max = 5;
    kappa = 1e16;
    GvacGGtemp = view(GvacGG, startidx1:startidx1+2, startidx2:startidx2+2);
    fill!(GvacGGtemp, zero(Complex{FT}));
    GvacGGSR = zeros(Complex{FT}, 3, 3);
    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));
    if (sigmatot < sqrt(Auc))
        kappa = sqrt(pi/Auc);
        N1max = 5;
        N2max = 5;
        setindex!(GvacGGSR, GFVACGGEwaldSR(freq, blochk, posfield, possource,
                                           latticevecs, alphafield, alphasource,
                                           kappa, N1max, N2max), 1:3, 1:3);
    end
    if (norm(posfield .- possource)/sigmatot <= eps())
        GvacGGtemp .+= GvacGGSR .+
        GFVACGGEwaldLR2D(freq, blochk, posfield, possource, reciprocalvecs, alphafield,
                         alphasource, kappa, M1max, M2max) -
        GFVACGGcoincident(freq, alphafield, alphasource)*I;
    else
        GvacGGtemp .+= GvacGGSR .+
        GFVACGGEwaldLR2D(freq, blochk, posfield, possource, reciprocalvecs, alphafield,
                         alphasource, kappa, M1max, M2max);
    end
    
end

"""

    GFVACGGEwaldSR(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                   posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                   latticevecs::AbstractArray{FT, 2},
                   alphafield::FT, alphasource::FT, kappa::FT,
                   N1max::Integer, N2max::Integer=0) where FT<:AbstractFloat

Return the short-range Ewald summation contribution to the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` and wavevector `blochk` (a
3-element vector), corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively, in 1 or 2 periodic
dimensions (with lattice parameters `latticevecs` for the real lattice vectors as
appropriate). The sums over the real lattice run over `-1*N1max:N1max` and `-1*N2max:N2max`
along each lattice vector, and `kappa` is the Ewald cutoff parameter.

"""
function GFVACGGEwaldSR(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                        posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                        latticevecs::AbstractArray{FT, 2},
                        alphafield::FT, alphasource::FT, kappa::FT,
                        N1max::Integer, N2max::Integer=0) where FT<:AbstractFloat

    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));
    
    q = sigmatot*freq/c;
    mynu = 1/hypot(1, 1/(sigmatot*kappa));
    
    n1list = -1*N1max:N1max;
    n2list = -1*N2max:N2max;
    
    GvacGGSR = zeros(Complex{FT}, 3, 3);
    currGvacGGSR = zeros(Complex{FT}, 3, 3);
    
    for n2 in n2list, n1 in n1list
      
        currlatticevec = n1*latticevecs[:,1] .+ n2*latticevecs[:,2];
        kdotR = dot(blochk, currlatticevec);
        myrhat = posfield .+ currlatticevec .- possource;
        rho = norm(myrhat)/sigmatot;
        nurho = mynu*rho;
        
        # simple expression for atoms coinciding in the same unit cell
        if (rho <= eps())
      
            setindex!(currGvacGGSR,
                      Matrix{eltype(currGvacGGSR)}((((2/sqrt(pi))*(q^2 - 1 - mynu*(q^2 - mynu^2)*exp((q^2 / 4) * (1/mynu^2 - 1))) + 1im*q^3 * (erf(1im*q/2) - erf(1im*q/(2*mynu))))/(6*pi*sigmatot^3))*I, 3, 3),
                      1:3, 1:3);
      
    # more complicated otherwise
        else
          
            fill!(currGvacGGSR, zero(eltype(currGvacGGSR)));
            normalize!(myrhat);
            term1 = exp(-1*(q^2 / 4) + 1im*q*rho)*erfc(-1im*(q/2) - rho);#exp(-0.25*q^2 + 1im*q*rho - (-0.5im*q - rho)^2 + log(erfcx(-0.5im*q - rho)));
            if (isnan(term1)) # just in case
                term1 = zero(eltype(currGvacGGSR));
            end
            term2 = exp(-1*(q^2 / 4) + 1im*q*rho)*erfc(-1im*q/(2*mynu) - nurho);#exp(-0.25*q^2 + 1im*q*rho - (-0.5im*q/mynu - nurho)^2 + log(erfcx(-0.5im*q/mynu - nurho)));
            if (isnan(term2)) # just in case
                term2 = zero(eltype(currGvacGGSR));
            end
# exp(-0.25*q^2 - 1im*q*rho)*erfc(-0.5im*q + rho), computed with better numerical stability given rho >= 0
            term3 = exp(-1*(q^2 / 4) - 1im*q*rho - (-1im*(q/2) + rho)^2 +
                        log(erfcx(-1im*(q/2) + rho)));
            if (isnan(term3)) # just in case
                term3 = zero(eltype(currGvacGGSR));
            end
# exp(-0.25*q^2 - 1im*q*rho)*erfc(-0.5im*q/mynu + nurho), computed with better numerical stability given rho >= 0
            term4 = exp(-1*(q^2 / 4) - 1im*q*rho - (-1im*q/(2*mynu) + nurho)^2 +
                        log(erfcx(-1im*q/(2*mynu) + nurho)));
            if (isnan(term4)) # just in case
                term4 = zero(eltype(currGvacGGSR));
            end

            gvacGGSR = ((term1 - term2) - (term3 - term4))/(8*pi*sigmatot*rho);
            if (isnan(gvacGGSR)) # just in case
              gvacGGSR = zero(gvacGGSR);
            end
            dgvacGGSRdrho = ((4/sqrt(pi))*(exp(-1*rho^2) * rho -
                                           exp((q^2 / 4) * (1/mynu^2 - 1) - (nurho)^2) *
                                           nurho) +
                             (-1 + 1im*q*rho)*(term1 - term2) +
                             (1 + 1im*q*rho)*(term3 - term4))/(8*pi*sigmatot*rho^2);
            if (isnan(dgvacGGSRdrho)) # just in case
                dgvacGGSRdrho = zero(dgvacGGSRdrho);
            end
            d2gvacGGSRdrho2 = ((8/sqrt(pi))*(exp((q^2 / 4) * (1/mynu^2 - 1) -
                                                 nurho^2) * (nurho + nurho^3) -
                                             exp(-1*rho^2) * (rho + rho^3)) -
                               (-2 + 2im*q*rho + (q*rho)^2)*(term1 - term2) +
                               (-2 - 2im*q*rho +
                                (q*rho)^2)*(term3 - term4))/(8*pi*sigmatot*rho^3);
            if (isnan(d2gvacGGSRdrho2)) # just in case
                d2gvacGGSRdrho2 = zero(d2gvacGGSRdrho2);
            end

            currGvacGGSR .+= exp(-1im*kdotR) .*
            (((d2gvacGGSRdrho2/(sigmatot^2) - dgvacGGSRdrho/(sigmatot^2 * rho))) .*
             complex(myrhat * transpose(myrhat)) .+
             Matrix{eltype(currGvacGGSR)}(((dgvacGGSRdrho/(sigmatot^2 * rho) +
                                            (freq/c)^2 * gvacGGSR))*I, 3, 3));
      
            if (sum(isnan.(currGvacGGSR)) > 0) # just in case
                fill!(currGvacGGSR, zero(eltype(currGvacGGSR)));
            end
        end

        GvacGGSR .+= currGvacGGSR;
        
    end
    
    return GvacGGSR;
    
end

"""

    GFVACGGEwaldLR1D(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                     posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                     reciprocalvec::AbstractArray{FT, 1},
                     alphafield::FT, alphasource::FT, kappa::FT,
                     Mmax::Integer, Smax::Integer) where FT<:AbstractFloat

Return the long-range Ewald summation contribution to the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` and wavevector `blochk` (a
3-element vector), corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively, in 1 periodic
dimension (with lattice parameters `reciprocalvec` for the reciprocal lattice vector as
appropriate). The sum over the reciprocal lattice runs over `-1*Mmax:Mmax`, the sum
computing the series representation of the reciprocal lattice expression runs over
`0:Smax`, and `kappa` is the Ewald cutoff parameter.

"""
function GFVACGGEwaldLR1D(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                          posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                          reciprocalvec::AbstractArray{FT, 1},
                          alphafield::FT, alphasource::FT, kappa::FT,
                          Mmax::Integer, Smax::Integer) where FT<:AbstractFloat

    if (real(freq) != 0.0 && imag(freq) != 0.0)
        error(string("GFVACGGEwaldLR1D(): freq is", freq,
                     "; GFVACGG for 1D-periodic systems can only handle freq purely",
                     " along the real or imaginary axes"));
    end

    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));

    q = sigmatot*freq/c;

    mynu2 = 1/(1 + 1/(sigmatot*kappa)^2);
    GvacGGLR = zeros(Complex{FT}, 3, 3);

# extract component of atomic displacement radial (perpendicular) with respect to array axis
    myreciprocalunitvec = normalize(reciprocalvec);
    myrvec = posfield .- possource;
    projperpreciprocal = myreciprocalunitvec * transpose(myreciprocalunitvec);
    subtractfromdiag!(projperpreciprocal, ones(FT, 3));
    myrhoperpvec = projperpreciprocal * myrvec/sigmatot;
    myrhoperp = norm(myrhoperpvec);
    nu2rhoperp2 = mynu2*myrhoperp^2;
    
    for mm=-1*Mmax:Mmax
      
        kgvec = blochk .+ reciprocalvec.*mm;
        h2over4nu2 = real(sigmatot^2 * dot(kgvec, kgvec) - q^2)/(4*mynu2);

        freqckg = Matrix{Complex{FT}}((freq/c)^2 * I, 3, 3) .-
        complex(kgvec * transpose(kgvec));

# if the source and field points lie along the same line, evaluate a simpler expression
        if (myrhoperp <= eps())
          
            currE1 = E1func(h2over4nu2);
            currE2 = exp(-1*h2over4nu2) - h2over4nu2*currE1;
            GvacGGLR .+= exp(-1*(q^2 / 4) + 1im*dot(kgvec, myrvec)) .*
            (freqckg .* currE1 .- (2*mynu2/sigmatot^2).*projperpreciprocal.*currE2);
            
        else # do the full sum over exponential integral orders

            sum1 = zero(FT);
            sum2 = zero(FT);
            sum3 = zero(FT);
            prevEint = E1func(h2over4nu2);
            currEint = copy(prevEint);
            coeff1 = zero(FT);
            coeff2 = zero(FT);
            coeff3 = zero(FT);
            
            for ss=0:Smax
              
                if (ss == 0)
                    coeff1 = one(coeff1);
                else
                
                    coeff1 = (-1/ss)*nu2rhoperp2 * coeff1;
                    currEint = (1/ss)*(exp(-1*h2over4nu2) - h2over4nu2*prevEint);
                    prevEint = copy(currEint);
                  
                    if (ss == 1)
                        coeff2 = -2.0*mynu2;
                    else
                    
                        coeff2 = (-1/(ss-1))*nu2rhoperp2 * coeff2;
                        if (ss == 2)
                            coeff3 = 4*mynu2^2;
                        else
                            coeff3 = (-1/(ss-2))*nu2rhoperp2 * coeff3;
                        end
                        
                    end
                    
                end
              
                sum1 += coeff1*currEint;
                sum2 += coeff2*currEint;
                sum3 += coeff3*currEint;
                
            end
            
            GvacGGLR .+= exp(-1*(q^2 / 4) + 1im*dot(kgvec, myrvec)) .*
            (sum1.*freqckg .+ sum2.*(projperpreciprocal./(sigmatot^2) .+
                                     (2im/sigmatot) .*
                                     sympart(complex(kgvec*transpose(myrhoperpvec)))) .+
             sum3.*complex(myrhoperpvec * transpose(myrhoperpvec))./(sigmatot^2));
            
        end
        
    end
    
    lmul!(norm(reciprocalvec)/(8*pi^2), GvacGGLR);
    return GvacGGLR;
    
end

"""

    GFVACGGEwaldLR2D(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                     posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                     reciprocalvecs::AbstractArray{FT, 2},
                     alphafield::FT, alphasource::FT, kappa::FT,
                     M1max::Integer, M2max::Integer) where FT<:AbstractFloat

Return the long-range Ewald summation contribution to the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq` and wavevector `blochk` (a
3-element vector), corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively, in 2 periodic
dimensions (with lattice parameters `reciprocalvecs` for the reciprocal lattice vectors as
appropriate). The sum over the reciprocal lattice runs over `-1*M1max:M1max` and
`-1*M2max:M2max` along each reciprocal lattice vector, and `kappa` is the Ewald cutoff
parameter.

"""
function GFVACGGEwaldLR2D(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                          posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                          reciprocalvecs::AbstractArray{FT, 2},
                          alphafield::FT, alphasource::FT, kappa::FT,
                          M1max::Integer, M2max::Integer) where FT<:AbstractFloat

    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));

    q = sigmatot*freq/c;

    mymu = sigmatot*kappa;
    mynu = mymu/hypot(1, mymu);
    mynu2 = 1/(1 + 1/mymu^2);
    GvacGGLR = zeros(Complex{FT}, 3, 3);
    
    ## unit vector normal to plane, and prefactor involving Brillouin zone area
    perpunitvec = cross(reciprocalvecs[:,1], reciprocalvecs[:,2]);
    ABZover16pi2 = norm(perpunitvec) / (16*pi^2);
    normalize!(perpunitvec);
    
    ## component of atomic displacement normal to the plane
    myrvec = posfield .- possource;
    myrhoperpsigned = dot(perpunitvec, myrvec)/sigmatot;
    
    for m2=-1*M2max:M2max, m1=-1*M1max:M1max
      
        kgvec = blochk .+ m1.*reciprocalvecs[:,1] .+ m2.*reciprocalvecs[:,2];
        kgnorm2 = dot(kgvec, kgvec);
        kgdotr = dot(kgvec, myrvec);
        
        freqckg = Matrix{Complex{FT}}((freq/c)^2 * I, 3, 3) .-
        complex(kgvec * transpose(kgvec));

        h = sqrt(complex(sigmatot^2 * kgnorm2 - q^2));
        
        ## if the source and field points do lie parallel to the lattice plane
        f = 2*exp(-1*(q^2 / 4) + 1im*kgdotr)*erfc(h/(2*mynu));
        dfdrhoperp = zero(f);
        d2fdrhoperp2 = h^2 * f -
        (4*mynu*h/sqrt(pi))*exp((q^2 / (4*mymu^2)) +
                                1im*kgdotr - sigmatot^2 * kgnorm2 / (4*mynu^2));
        
# if the source and field points do not lie parallel to the lattice plane
        if (abs(myrhoperpsigned) > eps())
          
# exp(-0.25*q^2 + 1im*kgdotr + h*myrhoperpsigned) * erfc(h/(2*mynu) + mynu*myrhoperpsigned), computed with possibly better numerical stability
            term1 = exp(-1*(q^2 / 4) +
                        1im*kgdotr + h*myrhoperpsigned -
                        (h/(2*mynu) + mynu*myrhoperpsigned)^2 +
                        log(erfcx(h/(2*mynu) + mynu*myrhoperpsigned)));
# no need for that in exponential decay
            term2 = exp(-1*(q^2 / 4) + 1im*kgdotr - h*myrhoperpsigned) *
            erfc(h/(2*mynu) - mynu*myrhoperpsigned);
#      f = exp(-0.25*q^2 + 1im*kgdotr + h*myrhoperpsigned) * erfc(h/(2*mynu) + mynu*myrhoperpsigned) + exp(-0.25*q^2 + 1im*kgdotr - h*myrhoperpsigned) * erfc(h/(2*mynu) - mynu*myrhoperpsigned);
            f = term1 + term2;
#      dfdrhoperp = h*(exp(-0.25*q^2 + 1im*kgdotr + h*myrhoperpsigned) * erfc(h/(2*mynu) - mynu*myrhoperpsigned) + exp(-0.25*q^2 + 1im*kgdotr - h*myrhoperpsigned) * erfc(h/(2*mynu) - mynu*myrhoperpsigned));
            dfdrhoperp = h*(term1 - term2);
            
            d2fdrhoperp2 = h^2 * f -
            (4*mynu*h/sqrt(pi))*
            exp((q^2 / (4*mymu^2)) + 1im*kgdotr -
                sigmatot^2 * kgnorm2/(4*mynu^2) - (mynu*myrhoperpsigned)^2);

        end
        
        GvacGGLR .+= (1/h) .* (sigmatot .* freqckg .* f .+
                               2im .* dfdrhoperp .*
                               sympart(complex(kgvec * transpose(perpunitvec))) .+
                               (1/sigmatot) .* d2fdrhoperp2 .*
                               complex(perpunitvec * transpose(perpunitvec)));

    end
    
    lmul!(ABZover16pi2, GvacGGLR);
    return GvacGGLR;
    
end

## The second set of functions in this file consists of aliases for
##  functions with fewer arguments: in particular, GFVACGG!(GvacGG,
##  startidx1, startidx2, freq, posfield, possource, alphafield,
##  alphasource) is the same as GFVACGG0!(GvacGG, startidx1,
##  startidx2, freq, posfield, possource, alphafield, alphasource),
##  which is convenient for compact molecules

"""

    GFVACGG!(GvacGG::AbstractArray{Complex{FT}, 2}, startidx1::Integer, startidx2::Integer,
             freq::Union{FT, Complex{FT}},
             posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
             alphafield::FT, alphasource::FT) where {FT<:AbstractFloat}

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with the 3-by-3 Green's
function interaction tensor in vacuum at frequency `freq`, corresponding to source position
`possource` and field position `posfield` (both 3-element vectors) described by Gaussian
basis functions parameterized by real atomic polarizabilities `alphasource` and `alphafield`
respectively, for nonperiodic geometries.

"""
GFVACGG!(GvacGG::AbstractArray{Complex{FT}, 2}, startidx1::Integer, startidx2::Integer,
         freq::Union{FT, Complex{FT}},
         posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
         alphafield::FT, alphasource::FT) where {FT<:AbstractFloat} =
GFVACGG0!(GvacGG, startidx1, startidx2, freq, posfield, possource, alphafield, alphasource);


"""

    GFVACGG_sca!(GvacGG::AbstractArray{Complex{FT}, 2}, startidx1::Integer, startidx2::Integer,
                 freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
                 possource::AbstractArray{FT, 1}, alphafield::FT,
                 alphasource::FT) where {FT<:AbstractFloat}

Fill `GvacGG[startidx1:startidx1+2, startidx2:startidx2+2]` in-place with a 3-by-3 block of
zeros, which is the scattering Green's function in vacuum by definition, irrespective of the
other function arguments. 

"""
GFVACGG_sca!(GvacGG::AbstractArray{Complex{FT}, 2}, startidx1::Integer, startidx2::Integer,
             freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
             possource::AbstractArray{FT, 1}, alphafield::FT,
             alphasource::FT) where {FT<:AbstractFloat} =
GFVACGG0_sca!(GvacGG, startidx1, startidx2, freq, posfield, possource, alphafield,
              alphasource);

## The third set of functions in this file consists of aliases for not
##  in-place function versions of the first and second sets


"""

    GFVACGG(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
            posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
            latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
            alphafield::FT, alphasource::FT, numdims::Integer) where FT<:AbstractFloat

Return the 3-by-3 Green's function interaction tensor in vacuum at frequency `freq` and
wavevector `blochk` (a 3-element vector), corresponding to source position `possource` and
field position `posfield` (both 3-element vectors) described by Gaussian basis functions
parameterized by real atomic polarizabilities `alphasource` and `alphafield` respectively,
using the appropriate function for `numdims` periodic dimensions (with lattice parameters
`latticevecs` for the real lattice vectors and `reciprocalvecs` for the reciprocal lattice
vectors as appropriate). `numdims` must be 0, 1, or 2.

"""
function GFVACGG(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                 posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                 latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
                 alphafield::FT, alphasource::FT, numdims::Integer) where FT<:AbstractFloat
  
    GvacGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFVACGG!(GvacGG, 1, 1, freq, blochk, posfield, possource,
             latticevecs, reciprocalvecs, alphafield, alphasource, numdims);
    return GvacGG;

end

"""

    GFVACGG_sca(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
                alphafield::FT, alphasource::FT,
                numdims::Integer) where {FT<:AbstractFloat}

Return a 3-by-3 block of zeros, which is the scattering Green's function in vacuum by
definition, irrespective of the other function arguments. 

"""
GFVACGG_sca(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
            posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
            latticevecs::AbstractArray{FT, 2}, reciprocalvecs::AbstractArray{FT, 2},
            alphafield::FT, alphasource::FT,
            numdims::Integer) where {FT<:AbstractFloat} = zeros(Complex{FT}, 3, 3);

"""

    GFVACGG0(freq::Union{FT, Complex{FT}},
             posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
             alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Return the 3-by-3 Green's function interaction tensor in vacuum at frequency `freq` without
spatial periodicity, corresponding to source position `possource` and field position
`posfield` (both 3-element vectors) described by Gaussian basis functions parameterized by
real atomic polarizabilities `alphasource` and `alphafield` respectively.

"""
function GFVACGG0(freq::Union{FT, Complex{FT}},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    GvacGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFVACGG0!(GvacGG, 1, 1, freq, posfield, possource, alphafield, alphasource);
    return GvacGG;
    
end

"""

    GFVACGG0_sca(freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
                 possource::AbstractArray{FT, 1}, alphafield::FT,
                 alphasource::FT) where {FT<:AbstractFloat}

Return a 3-by-3 block of zeros, which is the scattering Green's function in vacuum by
definition, irrespective of the other function arguments. 

"""
GFVACGG0_sca(freq::Union{FT, Complex{FT}}, posfield::AbstractArray{FT, 1},
             possource::AbstractArray{FT, 1}, alphafield::FT,
             alphasource::FT) where {FT<:AbstractFloat} = zeros(Complex{FT}, 3, 3);

"""

    GFVACGGcoincident(freq::Union{FT, Complex{FT}},
                      alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Return the scalar (as the tensor is the scalar multiplied by the 3-by-3 identity)
corresponding to vacuum Green's function interaction between two Gaussian basis functions,
defined by respective atomic polarizabilities `alphasource` and `alphafield`, when their
centers coincide; this quantity always has a finite imaginary part, and has a finite real
part when `alphafield` and `alphasource` are not both zero (meaning at least one basis
function has a finite width).

"""
function GFVACGGcoincident(freq::Union{FT, Complex{FT}},
                           alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    sigmatot = hypot(cbrt(alphafield), cbrt(alphasource)) / (cbrt(3) * sqrt(2*pi));
    q = sigmatot*freq/c;
    return (((2/sqrt(pi))*(q^2 - 1) + 1im*q^3 * erfc(-1im*q/2))/(6*pi*sigmatot^3));

end

"""

    GFVACGG1Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  latticevecs::AbstractArray{FT, 2},
                  reciprocalvecs::AbstractArray{FT, 2},
                  alphafield::FT, alphasource::FT) where FT<:AbstractFloat

Return the 3-by-3 Green's function interaction tensor in vacuum at frequency `freq` and
wavevector `blochk` (a 3-element vector), corresponding to source position `possource` and
field position `posfield` (both 3-element vectors) described by Gaussian basis functions
parameterized by real atomic polarizabilities `alphasource` and `alphafield` respectively,
in 1 periodic dimension (with lattice parameters `latticevecs` for the real lattice vectors
and `reciprocalvecs` for the reciprocal lattice vectors as appropriate). This is done via
Ewald summation over Gaussian basis functions in 1 dimension.

!!! warning

    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFVACGG1Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                       posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                       latticevecs::AbstractArray{FT, 2},
                       reciprocalvecs::AbstractArray{FT, 2},
                       alphafield::FT, alphasource::FT) where FT<:AbstractFloat

    GvacGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFVACGG1Ewald!(GvacGG, 1, 1, freq, blochk, posfield, possource,
                   latticevecs, reciprocalvecs, alphafield, alphasource);
    return GvacGG;

end

"""

    GFVACGG2Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                  posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                  latticevecs::AbstractArray{FT, 2},
                  reciprocalvecs::AbstractArray{FT, 2}, alphafield::FT,
                  alphasource::FT) where FT<:AbstractFloat

Return the 3-by-3 Green's function interaction tensor in vacuum at frequency `freq` and
wavevector `blochk` (a 3-element vector), corresponding to source position `possource` and
field position `posfield` (both 3-element vectors) described by Gaussian basis functions
parameterized by real atomic polarizabilities `alphasource` and `alphafield` respectively,
in 2 periodic dimensions (with lattice parameters `latticevecs` for the real lattice vectors
and `reciprocalvecs` for the reciprocal lattice vectors as appropriate). This is done via
Ewald summation over Gaussian basis functions in 2 dimensions.

!!! warning

    The Ewald summation convergence parameters are hard-coded into this function, though
    they should work for most common cases.

"""
function GFVACGG2Ewald(freq::Union{FT, Complex{FT}}, blochk::AbstractArray{FT, 1},
                       posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
                       latticevecs::AbstractArray{FT, 2},
                       reciprocalvecs::AbstractArray{FT, 2}, alphafield::FT,
                       alphasource::FT) where FT<:AbstractFloat

    GvacGG = Array{Complex{FT}, 2}(undef, 3, 3);
    GFVACGG2Ewald!(GvacGG, 1, 1, freq, blochk, posfield, possource,
                   latticevecs, reciprocalvecs, alphafield, alphasource);
    return GvacGG;

end

"""

    GFVACGG(freq::Union{FT, Complex{FT}},
            posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
            alphafield::FT, alphasource::FT) where {FT<:AbstractFloat}

Return the 3-by-3 Green's function interaction tensor in vacuum at frequency `freq`,
corresponding to source position `possource` and field position `posfield` (both 3-element
vectors) described by Gaussian basis functions parameterized by real atomic polarizabilities
`alphasource` and `alphafield` respectively, for nonperiodic geometries.

"""
GFVACGG(freq::Union{FT, Complex{FT}},
        posfield::AbstractArray{FT, 1}, possource::AbstractArray{FT, 1},
        alphafield::FT, alphasource::FT) where {FT<:AbstractFloat} =
GFVACGG0(freq, posfield, possource, alphafield, alphasource);
