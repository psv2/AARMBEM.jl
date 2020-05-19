"""

    constructMolAlpha!(myMolData::OneMol_WithPhonons{FT},
                       myPeriodicData::PeriodicData{FT},
                       freq::Union{FT, Complex{FT}},
                       k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the inverse molecular susceptibility `myMolData.alphainv` in-place at frequency
`freq` and wavevector `k`, the latter of which must be a 3-element real vector. With
phonons, `myMolData.alphainv` will be a matrix, the vector of atomic polarizabilities
`myMolData.alpha0` will depend on the elements of `myMolData.alphainv`, and
`myMolData.alphae` will be independent of these two. Without phonons, all three will be
simply related to each other.

"""
function constructMolAlpha!(myMolData::OneMol_WithPhonons{FT},
                            myPeriodicData::PeriodicData{FT},
                            freq::Union{FT, Complex{FT}},
                            k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat
  
    Ke3N = vecNto3N(myMolData.Ke);
    Qeinv3N = vecNto3N(1 ./ myMolData.Qe);
  
    setindex!(myMolData.alphainv, KIk(myMolData, myPeriodicData, k),
              1:3*myMolData.numatoms, 1:3*myMolData.numatoms);
    addtodiag!(myMolData.alphainv,
               Ke3N .- vecNto3N(1im.*freq.*myMolData.BI .+ freq^2 .* myMolData.MI));
    setindex!(myMolData.alphainv, inv(myMolData.alphainv),
              1:3*myMolData.numatoms, 1:3*myMolData.numatoms);
    mulMatDiagright!(myMolData.alphainv, Ke3N);
    mulMatDiagleft!(myMolData.alphainv, Ke3N);
    subtractfromdiag!(myMolData.alphainv,
                      Ke3N .- vecNto3N(1im.*freq.*myMolData.Be .+ freq^2 .* myMolData.Me));
    mulMatDiagright!(myMolData.alphainv, Qeinv3N);
    mulMatDiagleft!(myMolData.alphainv, Qeinv3N);

    alpha = inv(myMolData.alphainv);
    alphacarttr = Array{eltype(alpha), 2}(undef, myMolData.numatoms, myMolData.numatoms);
    @inbounds for qq=1:myMolData.numatoms
        @inbounds for pp=1:myMolData.numatoms
            alphacarttr[pp, qq] = tr(alpha[3*pp-2:3*pp, 3*qq-2:3*qq]);
        end
    end
  
    setindex!(myMolData.alpha0, abs.(vec(sum(alphacarttr, dims=2))), 1:myMolData.numatoms);
    lmul!(1/3, myMolData.alpha0);
    setindex!(myMolData.alphae,
              myMolData.Qe.^2 ./ (myMolData.Ke .-
                                  (1im.*freq.*myMolData.Be .+ freq^2 .* myMolData.Me)),
              1:myMolData.numatoms);
  
end

"""

    constructMolAlpha!(myMolData::OneMol_NoPhonons{FT},
                       myPeriodicData::PeriodicData{FT},
                       freq::Union{FT, Complex{FT}},
                       k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the inverse molecular susceptibility `myMolData.alphainv` in-place at frequency
`freq` and wavevector `k`, the latter of which must be a 3-element real vector. Without
phonons, `myMolData.alphainv` will be a vector, and the vectors `myMolData.alpha0` and
`myMolData.alphae` will be simply related to `myMolData.alphainv`.

"""
function constructMolAlpha!(myMolData::OneMol_NoPhonons{FT},
                            myPeriodicData::PeriodicData{FT},
                            freq::Union{FT, Complex{FT}},
                            k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    setindex!(myMolData.alphainv,
              (myMolData.Ke .- (1im.*freq.*myMolData.Be .+
                                freq^2 .* myMolData.Me))./(myMolData.Qe.^2),
              1:myMolData.numatoms);
    setindex!(myMolData.alphae, 1 ./ myMolData.alphainv, 1:myMolData.numatoms);
    setindex!(myMolData.alpha0, abs.(myMolData.alphae), 1:myMolData.numatoms);
  
end

"""

    KIk(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT},
        k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Construct the `k`-dependence of the internuclear spring constant matrix (only if phonons are
present), and return the result. Return `myMolData.KI` if the geometry has no periodic
dimensions.

See also: [`constructMolAlpha!(myMolData::OneMol_WithPhonons{FT},
                       myPeriodicData::PeriodicData{FT},
                       freq::Union{FT, Complex{FT}},
                       k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat`](@ref).

"""
function KIk(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT},
             k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    if (myPeriodicData.numdims < 1)
        return complex(myMolData.KI);
    end
    if (length(k) != 3)
        error("KIk(): k must have length 3");
    end
    currKIk = zeros(Complex{FT}, 3*myMolData.numatoms, 3*myMolData.numatoms);
    currKIblock = zeros(Complex{FT}, 3*myMolData.numatoms, 3*myMolData.numatoms);
    for nn2=1:myMolData.numblocks2
        for nn1=1:myMolData.numblocks1
            currKIblock[:, :] .= myMolData.KI[3*(nn1-1)*myMolData.numatoms+1:3*nn1*myMolData.numatoms,
                      3*(nn2-1)*myMolData.numatoms+1:3*nn2*myMolData.numatoms];
#            setindex!(currKIblock, myMolData.KI,
#                      3*(nn1-1)*myMolData.numatoms+1:3*nn1*myMolData.numatoms,
#                      3*(nn2-1)*myMolData.numatoms+1:3*nn2*myMolData.numatoms);
            lmul!(exp(-1im*dot(k,
                               myMolData.blocklist1[nn1] .*
                               myPeriodicData.latticevecs[:, 1] .+
                               myMolData.blocklist2[nn2] .*
                               myPeriodicData.latticevecs[:, 2])),
                  currKIblock);
            hermpart!(currKIblock);
            currKIk[:, :] .+= currKIblock[:, :];
        end
    end
    hermpart!(currKIk);
    return currKIk;
    
end

