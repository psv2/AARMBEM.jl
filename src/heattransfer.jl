"""

    heattransfer(myARGS)

Extract command-line arguments for computing the radiative heat transfer among a collection
of molecules, then compute the radiative heat transfer appropriately.

!!! note
    
    Any frequencies will be converted to their absolute values. Frequencies are not allowed
    to be exactly zero.

!!! note
    
    This computes the integrand
    ``\\operatorname{Tr}(\\operatorname{asym}(\\mathbb{V}_{m}^{-1}) \\mathbb{P}_{m} \\mathbb{T}^{\\dagger} \\operatorname{asym}(\\mathbb{P}_{n} \\mathbb{G}^{\\mathrm{env}}) \\mathbb{T}\\mathbb{P}_{m})``
    at specified real frequencies, without doing any integration. This includes neither
    numerical prefactors outside of the trace, nor the Planck function.

!!! tip

    This function can be embarrassingly parallelized over frequency and wavevector. However,
    the integrand at each frequency and wavevector is calculated for every transformation,
    as this maximizes the efficiency of the program by virtue of reusing quantities
    independent of geometric transformations.

"""
function heattransfer(myARGS)

    mollistfilename = readRequiredFilenameARG(myARGS, "mollistfilename", "AARMBEM_heat.jl");
    outfilename = readRequiredFilenameARG(myARGS, "outfilename", "AARMBEM_heat.jl");
  
    translistfilename = readOptionalFilenameARG(myARGS, "translistfilename");
    periodicfilename = readOptionalFilenameARG(myARGS, "periodicfilename");

    Genvstr = readGenvstr(myARGS);

    nophonons, outfilename = readOptionalBoolcondARG(myARGS, "nophonons", outfilename);
    nonretarded, outfilename = readOptionalBoolcondARG(myARGS, "nonretarded", outfilename);

    FT = readOptionalFloatType(myARGS);
    println("FT is ", FT);
    
    myAllMolData = readMolSystem(mollistfilename, nophonons, FT);
    myTransData = readTransData(myAllMolData, translistfilename);
    myPeriodicData = readPeriodicData(FT, periodicfilename);

    freqlist = Array{FT, 1}(undef, 0);
    readfreqlist!(freqlist, myARGS, "AARMBEM_heat.jl");

    klist = readklist(myARGS, "AARMBEM_heat.jl", myPeriodicData.numdims, FT);
    
    freqklist = hcat(vec(repeat(transpose(freqlist), size(klist, 1), 1)),
                     repeat(klist, length(freqlist), 1));

    heattransfer_general(outfilename, myAllMolData, myTransData, myPeriodicData, Genvstr,
                      nonretarded, freqklist)
  
end

"""

    heattransfer_general(outfilename, myAllMolData, myTransData, myPeriodicData, Genvstr,
                         nonretarded, freqklist)

Compute the actual radiative heat transfer integrand and print it to `outfilename`.

"""
function heattransfer_general(outfilename::AbstractString,
                              myAllMolData::MolSystem{<:OneMol{FT}},
                              myTransData::TransData{FT},
                              myPeriodicData::PeriodicData{FT},
                              Genvstr::AbstractString, nonretarded::Bool,
                              freqklist::Array{FT, 2}) where FT<:AbstractFloat

    GFSCAGG! = GFVACGG_sca!;
    if (Genvstr == "PEC")
        GFSCAGG! = GFPECGG!;
    end

    totalTinvmat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                            3*myAllMolData.cumnumatomslist[end]);
    totalGenvmat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                            3*myAllMolData.cumnumatomslist[end]);
    totalTmat = SharedArray{Complex{FT}}(3*myAllMolData.cumnumatomslist[end],
                                         3*myAllMolData.cumnumatomslist[end]);
    Phimn = zeros(FT, myAllMolData.nummols, myAllMolData.nummols);
    
    for ll=1:size(freqklist, 1)

        currfreq = abs(freqklist[ll, 1]);

        currfreqG = currfreq;
        if (nonretarded)
            currfreqG = zero(currfreqG);
        end
        currk = vec(freqklist[ll, 2:4]);
        fill!(totalTinvmat, zero(eltype(totalTinvmat)));
        fill!(totalGenvmat, zero(eltype(totalGenvmat)));
        fill!(totalTmat, zero(eltype(totalTmat)));

        for nn=1:myAllMolData.nummols

            currumolind = myAllMolData.uniquemolinds[nn];

            if (currumolind > -1)
            
                constructMolAlpha!(myAllMolData.uniqueMolDataArray[currumolind],
                                   myPeriodicData, currfreq, currk);
                
                assemble1MolGvacinf!(myAllMolData.uniqueMolDataArray[currumolind],
                                     myPeriodicData, currfreqG, currk);

                assemble1MolGvacinfdiag!(myAllMolData.uniqueMolDataArray[currumolind],
                                         currfreqG);
                
            end
            
        end


        for tt=1:myTransData.numTrans

            transformAtomPos!(myTransData, myAllMolData, tt);

            for nn2=1:myAllMolData.nummols
            
                startidx2 = 1;
                if (nn2 > 1)
                    startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
                end
                currmolind2 = myAllMolData.uniquemolinds[nn2];
                if (currmolind2 < 0)
                    currmolind2 = myAllMolData.dupemolinds[nn2];
                end

                if (myTransData.changedFromBefore[nn2, tt])
                    assembleTinvGenvDiagBlock!(totalTinvmat, totalGenvmat, myAllMolData,
                                               myTransData, nn2, currmolind2, startidx2,
                                               tt, myPeriodicData, GFSCAGG!, currfreqG,
                                               currk);
                end

            end

            for nn2=1:myAllMolData.nummols
            
                startidx2 = 1;
                if (nn2 > 1)
                    startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
                end

                for nn1=(nn2+1):myAllMolData.nummols

                    startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;

                    if (myTransData.changedFromBefore[nn2, tt] ||
                        myTransData.changedFromBefore[nn1, tt])

                        assembleTinvGenvOffDiagBlock!(totalTinvmat, totalGenvmat,
                                                      myAllMolData, myTransData, nn1, nn2,
                                                      startidx1, startidx2, myPeriodicData,
                                                      GFSCAGG!, currfreqG, currk);
                    end
                end

            end

            if (myPeriodicData.numdims > 0)
            
                for nn2=1:myAllMolData.nummols
                  
                    startidx2 = 1;
                    if (nn2 > 1)
                        startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
                    end

                    for nn1=(nn2+1):myAllMolData.nummols
                      
                        startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
                        
                        if (myTransData.changedFromBefore[nn2, tt] ||
                            myTransData.changedFromBefore[nn1, tt])
                          
                            assembleTinvGenvOffDiagBlock!(totalTinvmat, totalGenvmat,
                                                          myAllMolData, myTransData, nn2, nn1,
                                                          startidx2, startidx1, myPeriodicData,
                                                          GFSCAGG!, currfreqG, currk);
                        end
                    end

                end
            end
            
            totalTmat .= inv(totalTinvmat);
            
            assemblePhimn!(Phimn, myAllMolData, myPeriodicData, totalTmat, totalGenvmat);

            outfile = open(outfilename, "a+");
            if (myPeriodicData.numdims == 0)

                @printf(outfile, "%s %.16e", myTransData.transLabels[tt], abs(currfreq));
                for nn1=1:myAllMolData.nummols
                    for nn2=1:myAllMolData.nummols
                        @printf(outfile, " %.16e", Phimn[nn1, nn2]);;
                    end
                end
                @printf(outfile, "\n");
            else
                @printf(outfile, "%s %.16e %.16e %.16e %.16e", myTransData.transLabels[tt],
                        abs(currfreq), currk[1], currk[2], currk[3]);
                for nn1=1:myAllMolData.nummols
                    for nn2=1:myAllMolData.nummols
                        @printf(outfile, " %.16e", Phimn[nn1, nn2]);
                    end
                end
                @printf(outfile, "\n%s %.16e %.16e %.16e %.16e",
                        myTransData.transLabels[tt], abs(currfreq), -1*currk[1],
                        -1*currk[2], -1*currk[3]);
                for nn1=1:myAllMolData.nummols
                    for nn2=1:myAllMolData.nummols
                        @printf(outfile, " %.16e", Phimn[nn2, nn1]);
                    end
                end
                @printf(outfile, "\n");
            end
            close(outfile);
            
            untransformAtomPos!(myTransData, myAllMolData, tt);

        end
    
    end
  
end

"""

    assemblePhimn!(Phimn::Array{FT, 2},
                   myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                   myPeriodicData::PeriodicData{FT},
                   totalTmat::SharedArray{Complex{FT}, 2},
                   totalGenvmat::SharedArray{Complex{FT}, 2})

Assemble the matrix `Phimn`, whose elements `Phimn[nn1, nn2]` represent the energy transfer
from molecule `nn1` to molecule `nn2`, computed from trace expressions involving the
matrices `totalTmat` and `totalGenvmat`, in-place. If there is no spatial periodicity,
`Phimn` should be symmetric. If there is spatial periodicity, the transpose of `Phimn` at
wavevector ``\\vec{k}`` is the equivalent of computing `Phimn` at ``-\\vec{k}``. This assumes
the presence of phonons.

"""
function assemblePhimn!(Phimn::Array{FT, 2},
                        myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                        myPeriodicData::PeriodicData{FT},
                        totalTmat::SharedArray{Complex{FT}, 2},
                        totalGenvmat::SharedArray{Complex{FT}, 2}) where FT<:AbstractFloat

    
    for nn2=1:myAllMolData.nummols

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        Genvview = view(totalGenvmat,
                        (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                        1:3*myAllMolData.cumnumatomslist[end]);

        startmol1 = 1;
        if (myPeriodicData.numdims == 0)
            startmol1 = nn2;
        end
        
        for nn1=startmol1:myAllMolData.nummols
          
            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            Tview = view(totalTmat,
                         1:3*myAllMolData.cumnumatomslist[end],
                         (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]));
            
            currmolind1 = 1;
            if (myAllMolData.dupemolinds[nn1] == -1)
                currmolind1 = myAllMolData.uniquemolinds[nn1];
            else
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end
            Phimn[nn1, nn2] =
            imag(tr(ahermpart(myAllMolData.uniqueMolDataArray[currmolind1].alphainv) *
                    adjoint(Tview)[:, (startidx2-1) .+
                                   (1:3*myAllMolData.numatomslist[nn2])] *
                    Genvview * Tview));
            if (myPeriodicData.numdims == 0)
                Phimn[nn2, nn1] = Phimn[nn1, nn2];
            end
        end
        
    end

end
    
"""

    assemblePhimn!(Phimn::Array{FT, 2},
                   myAllMolData::MolSystem{OneMol_NoPhonons{FT}},
                   myPeriodicData::PeriodicData{FT},
                   totalTmat::SharedArray{Complex{FT}, 2},
                   totalGenvmat::SharedArray{Complex{FT}, 2})

Assemble the matrix `Phimn`, whose elements `Phimn[nn1, nn2]` represent the energy transfer
from molecule `nn1` to molecule `nn2`, computed from trace expressions involving the
matrices `totalTmat` and `totalGenvmat`, in-place. If there is no spatial periodicity,
`Phimn` should be symmetric. If there is spatial periodicity, the transpose of `Phimn` at
wavevector ``\\vec{k}`` is the equivalent of computing `Phimn` at ``-\\vec{k}``. This assumes
the absence of phonons.

"""
function assemblePhimn!(Phimn::Array{FT, 2},
                        myAllMolData::MolSystem{OneMol_NoPhonons{FT}},
                        myPeriodicData::PeriodicData{FT},
                        totalTinvmat::SharedArray{Complex{FT}, 2},
                        totalGenvmat::SharedArray{Complex{FT}, 2}) where FT<:AbstractFloat

    
    for nn2=1:myAllMolData.nummols

        startidx2 = 1;
        if (nn2 > 1)
            startidx2 = 3*myAllMolData.cumnumatomslist[nn2-1] + 1;
        end
        Genvview = view(totalGenvmat,
                        (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                        1:3*myAllMolData.cumnumatomslist[end]);

        startmol1 = 1;
        if (myPeriodicData.numdims == 0)
            startmol1 = nn2;
        end
        
        for nn1=startmol1:myAllMolData.nummols
          
            startidx1 = 1;
            if (nn1 > 1)
                startidx1 = 3*myAllMolData.cumnumatomslist[nn1-1] + 1;
            end
            Tview = view(totalTmat,
                         1:3*myAllMolData.cumnumatomslist[end],
                         (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]));
            
            currmolind1 = 1;
            if (myAllMolData.dupemolinds[nn1] == -1)
                currmolind1 = myAllMolData.uniquemolinds[nn1];
            else
                currmolind1 = myAllMolData.dupemolinds[nn1];
            end
            Phimn[nn1, nn2] =
            imag(tr(mulMatDiagleft(adjoint(Tview)[:, (startidx2-1) .+
                                                  (1:3*myAllMolData.numatomslist[nn2])] *
                                   Genvview * Tview,
                                   vecNto3N(imag.(myAllMolData.uniqueMolDataArray[currmolind1].alphainv)))));
            if (myPeriodicData.numdims == 0)
                Phimn[nn2, nn1] = Phimn[nn1, nn2];
            end
        end
        
    end

end

"""

    assembleTinvGenvDiagBlock!(totalTinvmat::SharedArray{Complex{FT}, 2},
                               totalGenvmat::SharedArray{Complex{FT}, 2},
                               myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                               myTransData::TransData{FT}, nn2::Integer,
                               currmolind2::Integer, startidx2::Integer,
                               tt::Integer, myPeriodicData::PeriodicData{FT},
                               GFSCAGG!::Function,
                               currfreqG::Union{FT, Complex{FT}},
                               currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular diagonal block of the scattering Green's function interaction matrix
of molecule `nn2` (whose unique molecular index `currmolind2` is precomputed and passed as
an argument) at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, then add the
vacuum contribution of molecule `nn2`, performing any required matrix rotations
corresponding to transformation `tt` for molecule `nn2`. Stamp that added to coincident
contributions into
`totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`,
or that (without coincident contributions) subtracted from the inverse molecular
susceptibility into
`totalTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Finally, undo any rotation corresponding to transformation `tt` in the inverse T-operator of
molecule `nn2` for reuse with future transformations. This assumes the presence of phonons.

"""
function assembleTinvGenvDiagBlock!(totalTinvmat::SharedArray{Complex{FT}, 2},
                                    totalGenvmat::SharedArray{Complex{FT}, 2},
                                    myAllMolData::MolSystem{OneMol_WithPhonons{FT}},
                                    myTransData::TransData{FT}, nn2::Integer,
                                    currmolind2::Integer, startidx2::Integer,
                                    tt::Integer, myPeriodicData::PeriodicData{FT},
                                    GFSCAGG!::Function,
                                    currfreqG::Union{FT, Complex{FT}},
                                    currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalGenvmat, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));

    fill!(tempview, zero(eltype(tempview)));

    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn2, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].alphainv,
                 myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                 myTransData.baselineRot[3*nn2-2:3*nn2, :]);

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                 myTransData.baselineRot[3*nn2-2:3*nn2, :]);
    
    tempview .+= myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf[:, :];

    totalTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                 (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])] .=
    myAllMolData.uniqueMolDataArray[currmolind2].alphainv .- tempview;

    addtodiag!(tempview,
               vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinfdiag));
    
    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].alphainv,
                 transpose(myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                           myTransData.baselineRot[3*nn2-2:3*nn2, :]));

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 transpose(myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                           myTransData.baselineRot[3*nn2-2:3*nn2, :]));

end

"""

    assembleTinvGenvDiagBlock!(totalTinvmat::SharedArray{Complex{FT}, 2},
                               totalGenvmat::SharedArray{Complex{FT}, 2},
                               myAllMolData::MolSystem{OneMol_NoPhonons{FT}},
                               myTransData::TransData{FT}, nn2::Integer,
                               currmolind2::Integer,
                               startidx2::Integer,
                               tt::Integer, myPeriodicData::PeriodicData{FT},
                               GFSCAGG!::Function,
                               currfreqG::Union{FT, Complex{FT}},
                               currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular diagonal block of the scattering Green's function interaction matrix
of molecule `nn2` (whose unique molecular index `currmolind2` is precomputed and passed as
an argument) at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, then add the
vacuum contribution of molecule `nn2`, performing any required matrix rotations
corresponding to transformation `tt` for molecule `nn2`. Stamp that added to coincident
contributions into
`totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`,
or that (without coincident contributions) subtracted from the inverse molecular
susceptibility into
`totalTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
Finally, undo any rotation corresponding to transformation `tt` in the inverse T-operator of
molecule `nn2` for reuse with future transformations. This assumes the absence of phonons.

"""
function assembleTinvGenvDiagBlock!(totalTinvmat::SharedArray{Complex{FT}, 2},
                                    totalGenvmat::SharedArray{Complex{FT}, 2},
                                    myAllMolData::MolSystem{OneMol_NoPhonons{FT}},
                                    myTransData::TransData{FT}, nn2::Integer,
                                    currmolind2::Integer,
                                    startidx2::Integer,
                                    tt::Integer, myPeriodicData::PeriodicData{FT},
                                    GFSCAGG!::Function,
                                    currfreqG::Union{FT, Complex{FT}},
                                    currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalGenvmat, (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));

    fill!(tempview, zero(eltype(tempview)));

    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn2, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                 myTransData.baselineRot[3*nn2-2:3*nn2, :]);
    
    tempview .+= myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf[:, :];

    totalTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                 (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])] .=
    subtractfromdiag(tempview,
                     vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].alphainv));

    addtodiag!(tempview,
               vecNto3N(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinfdiag));
    
    rotate3Nmat!(myAllMolData.uniqueMolDataArray[currmolind2].GvacGGinf,
                 transpose(myTransData.rotArray[3*nn2-2:3*nn2, 3*tt-2:3*tt] *
                           myTransData.baselineRot[3*nn2-2:3*nn2, :]));

end

"""

    assembleTinvGenvOffDiagBlock!(totalTinvmat::SharedArray{Complex{FT}, 2},
                                  totalGenvTinvmat::SharedArray{Complex{FT}, 2},
                                  myAllMolData::MolSystem{<:OneMol{FT}},
                                  myTransData::TransData{FT},
                                  nn1::Integer, nn2::Integer,
                                  startidx1::Integer,
                                  startidx2::Integer,
                                  myPeriodicData::PeriodicData{FT},
                                  GFSCAGG!::Function,
                                  currfreqG::Union{FT, Complex{FT}},
                                  currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

Assemble the molecular off-diagonal block of the scattering Green's function interaction
matrix of molecules `nn2` as the source and `nn1` as the field (whose respective unique
molecular indices `currmolind1` and `currmolind2` are precomputed and passed as arguments)
at frequency `currfreqG` and wavevector `currk` using `GFSCAGG!`, multiply on
the left by the electronic polarizabilities, add the vacuum contribution, and stamp that
into
`totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`
and its negative into
`totalTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]), (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])]`.
If there is no spatial periodicity, do the same for the (unconjugated) transpose, swapping
`nn2` and `nn1` without performing further Green's function computations.

"""
function assembleTinvGenvOffDiagBlock!(totalTinvmat::SharedArray{Complex{FT}, 2},
                                       totalGenvmat::SharedArray{Complex{FT}, 2},
                                       myAllMolData::MolSystem{<:OneMol{FT}},
                                       myTransData::TransData{FT},
                                       nn1::Integer, nn2::Integer,
                                       startidx1::Integer,
                                       startidx2::Integer,
                                       myPeriodicData::PeriodicData{FT},
                                       GFSCAGG!::Function,
                                       currfreqG::Union{FT, Complex{FT}},
                                       currk::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat

    tempview = view(totalGenvmat, (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                    (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]));
    fill!(tempview, zero(eltype(tempview)));

    assembleGvacMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, currfreqG, currk);
    assembleGscaMolBlock!(tempview, myAllMolData, myTransData, nn1, nn2, 1, 1,
                          myPeriodicData, GFSCAGG!, currfreqG, currk);

    totalTinvmat[(startidx1-1).+(1:3*myAllMolData.numatomslist[nn1]),
                 (startidx2-1).+(1:3*myAllMolData.numatomslist[nn2])] .=
    -1 .* tempview;

    if (myPeriodicData.numdims == 0)
        totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                     (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])] .=
        transpose(deepcopy(tempview));
        totalTinvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                     (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])] .= -1 .*
        totalGenvmat[(startidx2-1).+(1:3*myAllMolData.numatomslist[nn2]),
                     (startidx1-1).+(1:3*myAllMolData.numatomslist[nn1])];
    end

end
