"""

    assemble1MolGvacinf!(mySingleMolData::OneMol{FT}, myPeriodicData::PeriodicData{FT},
                         freq::Union{FT, Complex{FT}},
                         blochk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the Green's function interaction matrix of a single molecule in vacuum,
`mySingleMolData.GvacGGinf`, in-place, at frequency `freq` and wavevector `blochk` using
additional information from `myPeriodicData` to speed up the computation. In particular,
if the periodic dimensionality is 0, then this matrix is symmetric (whether `freq` is
imaginary or not). Even if the periodic dimensionality is not 0, if `freq` is strictly
imaginary, then this matrix is Hermitian. Both of those cases allow for computing only the
lower triangle of 3-by-3 blocks. Otherwise, the full matrix must be assembled.

!!! warning

    This function depends on `mySingleMolData.alpha0` having been initialized to the
    desired atomic susceptibilities at the same `freq` and `blochk`, but does not check
    this.

"""
function assemble1MolGvacinf!(mySingleMolData::OneMol{FT}, myPeriodicData::PeriodicData{FT},
                              freq::Union{FT, Complex{FT}},
                              blochk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    fill!(mySingleMolData.GvacGGinf, zero(eltype(mySingleMolData.GvacGGinf)));
    Gtemp = Array{eltype(mySingleMolData.GvacGGinf), 2}(undef, 3, 3);
    if (myPeriodicData.numdims == 0)
      
        ## @inbounds @sync @distributed for qq=1:mySingleMolData.numatoms
        @inbounds for qq=1:mySingleMolData.numatoms
            @inbounds for pp=qq:mySingleMolData.numatoms

                GFVACGG!(Gtemp, 1, 1, freq, view(mySingleMolData.atompos0, 1:3, pp),
                        view(mySingleMolData.atompos0, 1:3, qq),
                        mySingleMolData.alpha0[pp], mySingleMolData.alpha0[qq]);
                mySingleMolData.GvacGGinf[3*pp-2:3*pp, 3*qq-2:3*qq] .= Gtemp;
                mySingleMolData.GvacGGinf[3*qq-2:3*qq, 3*pp-2:3*pp] .= transpose(Gtemp);

            end
        end

    elseif (real(freq) == zero(real(freq)))

        ## @inbounds @sync @distributed for qq=1:mySingleMolData.numatoms
        @inbounds for qq=1:mySingleMolData.numatoms
            @inbounds for pp=qq:mySingleMolData.numatoms
            
                GFVACGG!(Gtemp, 1, 1, freq, blochk, view(mySingleMolData.atompos0, 1:3, pp),
                         view(mySingleMolData.atompos0, 1:3, qq),
                         myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                         mySingleMolData.alpha0[pp], mySingleMolData.alpha0[qq],
                         myPeriodicData.numdims);
                mySingleMolData.GvacGGinf[3*pp-2:3*pp, 3*qq-2:3*qq] .= Gtemp;
                mySingleMolData.GvacGGinf[3*qq-2:3*qq, 3*pp-2:3*pp] .=
                adjoint(Gtemp);
            end
        end

    else

        ## @inbounds @sync @distributed for qq=1:mySingleMolData.numatoms
        @inbounds for qq=1:mySingleMolData.numatoms
            @inbounds for pp=1:mySingleMolData.numatoms

                GFVACGG!(Gtemp, 1, 1, freq, blochk,
                         view(mySingleMolData.atompos0, 1:3, pp),
                         view(mySingleMolData.atompos0, 1:3, qq),
                         myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                         mySingleMolData.alpha0[pp], mySingleMolData.alpha0[qq],
                         myPeriodicData.numdims);
                mySingleMolData.GvacGGinf[3*pp-2:3*pp, 3*qq-2:3*qq] .= Gtemp;
                
            end
        end

    end

end

"""

    assemble1MolGvacinfdiag!(mySingleMolData::OneMol{FT},
                             freq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

Assemble the scalar coincident Green's function interaction matrix elements of a single
molecule in vacuum, `mySingleMolData.GvacGGinfdiag`, in-place, at frequency `freq`.

!!! warning

    This function depends on `mySingleMolData.alpha0` having been initialized to the
    desired atomic susceptibilities at the same `freq` (and wavevector), but does not check
    this.

"""
function assemble1MolGvacinfdiag!(mySingleMolData::OneMol{FT},
                                  freq::Union{FT, Complex{FT}}) where FT<:AbstractFloat

    mySingleMolData.GvacGGinfdiag .= map(u -> GFVACGGcoincident(freq, u, u),
                                         mySingleMolData.alpha0);

end

"""

    subtract1MolalphainvGvac!(myMolData::OneMol{<:AbstractFloat})

Subtract the Green's function interaction matrix of a single molecule in vacuum,
`myMolData.GvacGGinf`, in-place, from its inverse susceptibility `myMolData.alphainv`.
Separately dispatch cases with versus without phonons.

!!! warning

    This function depends on `myMolData.alphainv` having been initialized to the desired
    molecular susceptibility at the same `freq` and `blochk`, but does not check this.

"""
function subtract1MolalphainvGvac!(myMolData::OneMol_NoPhonons{<:AbstractFloat})

    subtractfromdiag!(myMolData.GvacGGinf, vecNto3N(myMolData.alphainv));

end
function subtract1MolalphainvGvac!(myMolData::OneMol_WithPhonons{<:AbstractFloat})

    myMolData.GvacGGinf[:, :] .= myMolData.alphainv[:, :] .- myMolData.GvacGGinf[:, :];

end

"""

    assembleGscaMolBlock!(A::SharedArray{Complex{FT}, 2},
                          myAllMolData::MolSystem{<:OneMol{FT}},
                          myTransData::TransData{FT}, nn1::Integer, nn2::Integer,
                          startidx1::Integer, startidx2::Integer,
                          myPeriodicData::PeriodicData{FT},
                          GFSCAGG!::Function, freq::Union{FT, Complex{FT}},
                          blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the scattering Green's function interaction matrix (evaluated via the function
`GFSCAGG!`) between molecules `nn1` and `nn2` (sources evaluated in the latter, fields
evaluated in the former), at frequency `freq` and wavevector `blochk`, using atomic
polarizability and coordinate data respectively from `myAllMolData` and `myTransData`, and
add the result in-place to the matrix block
`A[startidx1:endidx1, startidx2:endidx2]` where `endidx1 - startidx1 + 1` is 3 times the
number of atoms in molecule `nn1`, and likewise `endidx2 - startidx2 + 1` is 3 times the
number of atoms in molecule `nn2`. If `nn1` and `nn2` are the same, the matrix
block is symmetric if the periodic dimensionality is 0, or Hermitian if the periodic
dimensionality is nonzero but `freq` is purely imaginary, so only the lower triangle of
blocks is assembled. Otherwise, the full set of matrix blocks is assembled. 

!!! note
    
    If `GFSCAGG!` is the scattering Green's function of vacuum (returning a 3-by-3 block of
    zeros), namely `GFVACGG_sca!` or `GFVACGG0_sca!`, this function does nothing, saving on
    computations.

!!! warning
    
    Currently, only `GFPECGG!` and `GFPECGG0!` are supported for `GFSCAGG!`, but this is not
    explicitly checked, and may change in future versions of this code.

!!! warning

    This function depends on `mySingleMolData.alpha0` having been initialized to the desired
    atomic susceptibilities at the same `freq` and `blochk`, as well as
    `myTransData.allAtomPos` having been modified to reflect the desired transformation, but
    does not check these.

"""
function assembleGscaMolBlock!(A::AbstractArray{Complex{FT}, 2},
                               myAllMolData::MolSystem{<:OneMol{FT}},
                               myTransData::TransData{FT}, nn1::Integer, nn2::Integer,
                               startidx1::Integer, startidx2::Integer,
                               myPeriodicData::PeriodicData{FT},
                               GFSCAGG!::Function, freq::Union{FT, Complex{FT}},
                               blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    if ((GFSCAGG!) != (GFVACGG_sca!) && (GFSCAGG!) != (GFVACGG0_sca!)) # only do anything if we aren't in vacuum
 
        currmolind1 = myAllMolData.dupemolinds[nn1];
        currmolind2 = myAllMolData.dupemolinds[nn2];
        if (currmolind1 < 0)
            currmolind1 = myAllMolData.uniquemolinds[nn1];
        end
        if (currmolind2 < 0)
            currmolind2 = myAllMolData.uniquemolinds[nn2];
        end

        if (size(A, 1) < 3*myAllMolData.numatomslist[nn1])
            error("A (currently with ", size(A, 1), " rows) should have at least ",
                  3*myAllMolData.numatomslist[nn1], " rows in assembleGenvMolBlock!()");
        end
        if (size(A, 2) < 3*myAllMolData.numatomslist[nn2])
            error("A (currently with ", size(A, 2), " columns) should have at least ",
                  3*myAllMolData.numatomslist[nn2], " columns in assembleGenvMolBlock!()");
        end
        if (startidx1 < 1 || startidx1 > size(A, 1) - 3*myAllMolData.numatomslist[nn1] + 1)
            error("assembleGenvMolBlock!(): startidx1 (currently ", startidx1,
                  ") must be between 1 and ",
                  size(A, 1) - 3*myAllMolData.numatomslist[nn1] + 1,
                  " as A has ", size(A, 1), " rows and the block will have ",
                  3*myAllMolData.numatomslist[nn1], " rows");
        end
        if (startidx2 < 1 || startidx2 > size(A, 2) - 3*myAllMolData.numatomslist[nn2] + 1)
            error("assembleGenvMolBlock!(): startidx2 (currently ", startidx2,
                  ") must be between 1 and ",
                  size(A, 2) - 3*myAllMolData.numatomslist[nn2] + 1,
                  " as A has ", size(A, 2), " columns and the block will have ",
                  3*myAllMolData.numatomslist[nn2], " columns");
        end
  
        Gtemp = Array{eltype(A), 2}(undef, 3, 3);
        fill!(Gtemp, zero(eltype(Gtemp)));
        allatomoffset1 = 0;
        if (nn1 > 1)
            allatomoffset1 = myAllMolData.cumnumatomslist[nn1-1];
        end
        allatomoffset2 = 0;
        if (nn2 > 1)
            allatomoffset2 = myAllMolData.cumnumatomslist[nn2-1];
        end
        posfield = zeros(eltype(myTransData.allAtomPos), 3);
        possource = zeros(eltype(myTransData.allAtomPos), 3);
        
        if (nn1 == nn2 && myPeriodicData.numdims == 0)

            @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
                possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
                posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + qq];
                GFSCAGG!(Gtemp, 1, 1, freq, posfield, possource,
                         myAllMolData.uniqueMolDataArray[currmolind1].alpha0[qq],
                         myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq]);
                sympart!(Gtemp);
                A[(startidx1-1).+(3*qq-2:3*qq), (startidx2-1).+(3*qq-2:3*qq)] .+=
                    Gtemp;
            end

            @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
                possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
                @inbounds for pp=(qq+1):myAllMolData.numatomslist[currmolind1]
                    posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + pp];
                    GFSCAGG!(Gtemp, 1, 1, freq, posfield, possource,
                             myAllMolData.uniqueMolDataArray[currmolind1].alpha0[pp],
                             myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq]);
                    A[(startidx1-1).+(3*pp-2:3*pp), (startidx2-1).+(3*qq-2:3*qq)] .+=
                    Gtemp;
                    A[(startidx2-1).+(3*qq-2:3*qq), (startidx1-1).+(3*pp-2:3*pp)] .+=
                    transpose(Gtemp);
                end
            end

        elseif (nn1 == nn2 && real(freq) == zero(real(freq)))

            @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
                possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
                posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + qq];
                GFSCAGG!(Gtemp, 1, 1, freq, blochk, posfield, possource,
                         myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                         myAllMolData.uniqueMolDataArray[currmolind1].alpha0[qq],
                         myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq],
                         myPeriodicData.numdims);
                hermpart!(Gtemp);
                A[(startidx1-1).+(3*qq-2:3*qq), (startidx2-1).+(3*qq-2:3*qq)] .+=
                    Gtemp;                        
            end

            @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
                possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
                @inbounds for pp=(qq+1):myAllMolData.numatomslist[currmolind1]
                    posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + pp];
                    GFSCAGG!(Gtemp, 1, 1, freq, blochk, posfield, possource,
                             myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                             myAllMolData.uniqueMolDataArray[currmolind1].alpha0[pp],
                             myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq],
                             myPeriodicData.numdims);
                    A[(startidx1-1).+(3*pp-2:3*pp), (startidx2-1).+(3*qq-2:3*qq)] .+=
                    Gtemp;
                    A[(startidx2-1).+(3*qq-2:3*qq), (startidx1-1).+(3*pp-2:3*pp)] .+=
                    adjoint(Gtemp);
                end
            end

        else

            @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
                possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
                @inbounds for pp=1:myAllMolData.numatomslist[currmolind1]
                    posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + pp];
                    GFSCAGG!(Gtemp, 1, 1, freq, blochk, posfield, possource,
                             myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                             myAllMolData.uniqueMolDataArray[currmolind1].alpha0[pp],
                             myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq],
                             myPeriodicData.numdims);
                    A[(startidx1-1).+(3*pp-2:3*pp), (startidx2-1).+(3*qq-2:3*qq)] .+=
                    Gtemp;
                end
            end

        end

    end
  
end
  
"""

    assembleGvacMolBlock!(A::AbstractArray{Complex{FT}, 2},
                          myAllMolData::MolSystem{<:OneMol{FT}},
                          myTransData::TransData{FT}, nn1::Integer, nn2::Integer,
                          startidx1::Integer, startidx2::Integer,
                          myPeriodicData::PeriodicData{FT},
                          freq::Union{FT, Complex{FT}},
                          blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the vacuum Green's function interaction matrix (evaluated via the function
`GFVACGG!`) between molecules `nn1` and `nn2` (sources evaluated in the latter, fields
evaluated in the former), at frequency `freq` and wavevector `blochk`, using atomic
polarizability and coordinate data respectively from `myAllMolData` and `myTransData`, and
add the result in-place to the matrix block
`A[startidx1:endidx1, startidx2:endidx2]` where `endidx1 - startidx1 + 1` is 3 times the
number of atoms in molecule `nn1`, and likewise `endidx2 - startidx2 + 1` is 3 times the
number of atoms in molecule `nn2`. If `nn1` and `nn2` are the same, the matrix
block is symmetric if the periodic dimensionality is 0, or Hermitian if the periodic
dimensionality is nonzero but `freq` is purely imaginary, so only the lower triangle of
blocks is assembled. Otherwise, the full set of matrix blocks is assembled. 

!!! tip

    Use `assemble1MolGvacinf!` when `nn1` and `nn2` are the same, then stamp that into `A`
    as needed. Only use this when `nn1` and `nn2` are different.

!!! warning

    This function depends on `mySingleMolData.alpha0` having been initialized to the desired
    atomic susceptibilities at the same `freq` and `blochk`, as well as
    `myTransData.allAtomPos` having been modified to reflect the desired transformation, but
    does not check these.

"""
function assembleGvacMolBlock!(A::AbstractArray{Complex{FT}, 2},
                               myAllMolData::MolSystem{<:OneMol{FT}},
                               myTransData::TransData{FT},
                               nn1::Integer, nn2::Integer,
                               startidx1::Integer,
                               startidx2::Integer,
                               myPeriodicData::PeriodicData{FT},
                               freq::Union{FT, Complex{FT}},
                               blochk::AbstractArray{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    currmolind1 = myAllMolData.dupemolinds[nn1];
    if (currmolind1 < 0)
        currmolind1 = myAllMolData.uniquemolinds[nn1];
    end
    currmolind2 = myAllMolData.dupemolinds[nn2];
    if (currmolind2 < 0)
        currmolind2 = myAllMolData.uniquemolinds[nn2];
    end

    if (size(A, 1) < 3*myAllMolData.numatomslist[nn1])
        error("A (currently with ", size(A, 1), " rows) should have at least ",
              3*myAllMolData.numatomslist[nn1], " rows in assembleGenvMolBlock!()");
    end
    if (size(A, 2) < 3*myAllMolData.numatomslist[nn2])
        error("A (currently with ", size(A, 2), " rows) should have at least ",
              3*myAllMolData.numatomslist[nn2], " rows in assembleGenvMolBlock!()");
    end
    if (startidx1 < 1 || startidx1 > size(A, 1) - 3*myAllMolData.numatomslist[nn1] + 1)
        error("assembleGenvMolBlock!(): startidx1 (currently ", startidx1,
              ") must be between 1 and ",
              size(A, 1) - 3*myAllMolData.numatomslist[nn1] + 1,
              " as A has ", size(A, 1), " rows and the block will have ",
              3*myAllMolData.numatomslist[nn1], " rows");
    end
    if (startidx2 < 1 || startidx2 > size(A, 2) - 3*myAllMolData.numatomslist[nn2] + 1)
        error("assembleGenvMolBlock!(): startidx2 (currently ", startidx2,
              ") must be between 1 and ",
              size(A, 2) - 3*myAllMolData.numatomslist[nn2] + 1,
              " as A has ", size(A, 2), " columns and the block will have ",
              3*myAllMolData.numatomslist[nn2], " columns");
    end

    Gtemp = Array{eltype(A), 2}(undef, 3, 3);
    fill!(Gtemp, zero(eltype(Gtemp)));
    allatomoffset1 = 0;
    if (nn1 > 1)
        allatomoffset1 = myAllMolData.cumnumatomslist[nn1-1];
    end
    allatomoffset2 = 0;
    if (nn2 > 1)
        allatomoffset2 = myAllMolData.cumnumatomslist[nn2-1];
    end
    posfield = zeros(eltype(myTransData.allAtomPos), 3);
    possource = zeros(eltype(myTransData.allAtomPos), 3);

    if (nn1 == nn2 && myPeriodicData.numdims == 0)
      
        @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
            possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
            @inbounds for pp=(qq+1):myAllMolData.numatomslist[currmolind1]
                posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + pp];
                GFVACGG!(Gtemp, 1, 1, freq, posfield, possource,
                         myAllMolData.uniqueMolDataArray[currmolind1].alpha0[pp],
                         myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq]);
                A[(startidx1-1).+(3*pp-2:3*pp), (startidx2-1).+(3*qq-2:3*qq)] .+=
                Gtemp;
                A[(startidx2-1).+(3*qq-2:3*qq), (startidx1-1).+(3*pp-2:3*pp)] .+=
                transpose(Gtemp);
            end
        end

    elseif (nn1 == nn2 && real(freq) == zero(real(freq)))

        @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
            possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
            @inbounds for pp=(qq+1):myAllMolData.numatomslist[currmolind1]
                posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + pp];
                GFVACGG!(Gtemp, 1, 1, freq, blochk, posfield, possource,
                         myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                         myAllMolData.uniqueMolDataArray[currmolind1].alpha0[pp],
                         myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq],
                         myPeriodicData.numdims);
                A[(startidx1-1).+(3*pp-2:3*pp), (startidx2-1).+(3*qq-2:3*qq)] .+=
                Gtemp;
                A[(startidx2-1).+(3*qq-2:3*qq), (startidx1-1).+(3*pp-2:3*pp)] .+=
                adjoint(Gtemp);
            end
        end

    else

        @inbounds for qq=1:myAllMolData.numatomslist[currmolind2]
            possource .= myTransData.allAtomPos[1:3, allatomoffset2 + qq];
            @inbounds for pp=1:myAllMolData.numatomslist[currmolind1]
                posfield .= myTransData.allAtomPos[1:3, allatomoffset1 + pp];
                GFVACGG!(Gtemp, 1, 1, freq, blochk, posfield, possource,
                         myPeriodicData.latticevecs, myPeriodicData.reciprocalvecs,
                         myAllMolData.uniqueMolDataArray[currmolind1].alpha0[pp],
                         myAllMolData.uniqueMolDataArray[currmolind2].alpha0[qq],
                         myPeriodicData.numdims);
                A[(startidx1-1).+(3*pp-2:3*pp), (startidx2-1).+(3*qq-2:3*qq)] .+=
                Gtemp;
            end
        end

    end

end
  
