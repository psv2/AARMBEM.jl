"""

    TransData{FT<:AbstractFloat}

Store information relevant to the particular geometrical configuration, not intrinsic to
molecules themselves, involving rigid geometrical transformations (translations and
rotations); this depends on the number of molecules `nummols` from
`MolSystem{T<:OneMol{<:AbstractFloat}}`.

`numTrans` is the number of transformations ``t`` (an integer).

`transArray` is a `3*nummols`-by-`numTrans` real array where each 3-element vector
`transArray[3*nn-2:3*nn, tt]` is the center-of-mass translation vector (relative to the
baseline) for molecule `nn` in transformation labeled `tt`.

`rotArray` is a `3*nummols`-by-`3*numTrans` real array where each 3-by-3 block
`rotArray[3*nn-2:3*nn, 3*tt-2:3*tt]` is the center-of-mass rotation matrix (relative to the
baseline) for molecule `nn` in transformation labeled `tt`.

`transLabels` is a `numTrans`-element vector of strings labeling each transformation `tt`.

`changedFromBefore` is an `nummols`-by-`numTrans` `BitArray` such that
`changedFromBefore[nn, tt]` is `true` if molecule `nn` is transformed in transformation `tt`
relative to transformation `tt - 1`, `false` otherwise.

`allAtomPos` is the horizontal concatenation of atomic coordinates (initially unmodified,
later modified, by transformations) across all molecules.

`baselineTrans` is a `3*nummols`-element real vector where each 3-element vector
`baselineTrans[3*nn-2:3*nn]` is the baseline center of mass coordinate for molecule `nn`.

`baselineRot` is a `3*nummols`-by-3 real matrix where each 3-by-3 block
`baselineRot[3*nn-2:3*nn,:]` is the baseline center of mass rotation matrix for the
orientation of molecule `nn`. Note that for both `baselineRot` and `rotArray`, each 3-by-3
Cartesian rotation matrix is given as the transpose of the overall matrix
``Z_{1} Y_{2} X_{3}`` in the Tait--Bryan convention, where angle 1 is ``\\psi``, angle 2 is
``\\theta``, and angle 3 is ``\\varphi``. In particular, each rotation matrix is the product
``R_{z}(-\\psi) R_{y}(-\\theta) R_{x}(-\\varphi)``, where ``R_{i}(\\alpha)`` is a standard
alibi rotation about Cartesian axis ``\\vec{e}_{i}``; for example, ``R_{z}(\\pi/2)``
multiplied on the left of the 3-element column ``(1, 0, 0)^{\\top}`` yields
``(0, 1, 0)^{\\top}``, and on the left of ``(0, 1, 0)^{\\top}`` yields
``(-1, 0, 0)^{\\top}``. Note also that in the file with name `translistfilename`, the angles
for each baseline orientation as well as each transformation beyond that must be listed in
the order `phi` then `theta` then `psi`, in DEGREES.
"""
struct TransData{FT<:AbstractFloat}

    transArray::Array{FT, 2}
    rotArray::Array{FT, 2}
    transLabels::Array{<:AbstractString, 1}
    changedFromBefore::BitArray{2}
    numTrans::Integer
    allAtomPos::Array{FT, 2}
    baselineTrans::Array{FT, 1}
    baselineRot::Array{FT, 2}
  
end

"""

    readTransData(myMolData::MolSystem{<:OneMol{FT}},
                  translistfilename::AbstractString="") where FT<:AbstractFloat

Construct TransData from existing molecular data `myMolData` and a file of name
`translistfilename` containing information about geometric transformations; if no file is
specified, produce a default set of transformations (no displacements or rotations at all).

See also: [`MolSystem{<:OneMol{FT}}`](@ref)

"""
function readTransData(myMolData::MolSystem{<:OneMol{FT}},
                       translistfilename::AbstractString="") where FT<:AbstractFloat

    ## default: no transformation file specified
    if (translistfilename == "")
        numtrans = 1;
        allmoltrans = zeros(FT, numtrans, 6*myMolData.nummols);
        allmoltranslabels = ["DEFAULT"];
        changedfrombefore = falses(myMolData.nummols, numtrans);
        transArray = zeros(FT, 3*myMolData.nummols, numtrans);
        rotArray = repeat(Matrix{FT}(I, 3, 3), myMolData.nummols, numtrans);
        
        allAtomPos = Array{FT,2}(undef, 3, myMolData.cumnumatomslist[end]);
        concatenateAtomPos!(myMolData, allAtomPos);
        
        baselineTrans = zeros(FT, 3*myMolData.nummols);
        baselineRot = repeat(Matrix{FT}(I, 3, 3), myMolData.nummols, 1);
        
        myTransData = TransData{FT}(transArray, rotArray, allmoltranslabels,
                                    changedfrombefore, numtrans, allAtomPos,
                                    baselineTrans, baselineRot);
        initializeAtomPos!(myTransData, myMolData);
        return myTransData;
    end

    # read file and format each line to get rid of extraneous spaces and things like that
    transstrarray = strip.(chomp.(readlines(translistfilename)));
    transstartinds = findall(map(u -> startswith(u, "TRANS "), transstrarray));
    transendinds = findall(transstrarray .== "ENDTRANS");
    numtrans = length(transendinds);
    if (length(transstartinds) != numtrans)
        error(string(translistfilename, ": number of TRANS and ENDTRANS statements",
                     " should be the same"));
    end
    transbothinds = vec(transpose(hcat(transstartinds, transendinds)));
    if (transbothinds != sort(transbothinds))
        error(string(translistfilename, ": TRANS and ENDTRANS statements should",
                     " alternate in ordered occurrence"));
    end

    # populate the lists of transformation labels, translation vectors, and Euler angles
    allmoltranslabels = Array{String, 1}(undef, numtrans);
    allmoltrans = zeros(FT, 6*myMolData.nummols + 1, numtrans);
    for tt=1:numtrans
        allmoltranslabels[tt] = replace(strip(chomp(transstrarray[transstartinds[tt]])),
                                        "TRANS " => "");
        currtransstrarray = strip.(chomp.(transstrarray[transstartinds[tt]+1:transendinds[tt]-1]));
        for currtransstr in currtransstrarray
            currtransstrsplit = split(currtransstr);
            currmollabel = currtransstrsplit[1];
            currmolind = findfirst(myMolData.molLabels .== currmollabel);

            if (currmolind != nothing)
                ## assemble transformation list in the following order for each row in a transformation:
                ## [phi_1, theta_1, psi_1, ..., phi_N, theta_N, psi_N, x_1, y_1, z_1, ..., x_N, y_N, z_N]
                ## Note: all Euler angles should be given in the transformation file as ALIBI transformations in DEGREES
                allmoltrans[3*myMolData.nummols+3*currmolind-2,
                            tt] = parse(FT,currtransstrsplit[2]);
                allmoltrans[3*myMolData.nummols+3*currmolind-1,
                            tt] = parse(FT, currtransstrsplit[3]);
                allmoltrans[3*myMolData.nummols+3*currmolind,
                            tt] = parse(FT, currtransstrsplit[4]);
                if (length(currtransstrsplit) == 7)
                    allmoltrans[3*currmolind-2, tt] = parse(FT, currtransstrsplit[5]);
                    allmoltrans[3*currmolind-1, tt] = parse(FT, currtransstrsplit[6]);
                    allmoltrans[3*currmolind, tt] = parse(FT, currtransstrsplit[7]);
                end
            end
        end
        allmoltrans[end, tt] = tt;
    end

    # populate the lists of baseline translation vectors and Euler angles
    baselinestartind = findfirst(map(u -> startswith(u, "BASELINE"), transstrarray));
    baselineendind = findfirst(map(u -> startswith(u, "ENDBASELINE"), transstrarray));
    baselinetrans = zeros(FT, 6*myMolData.nummols);

    if (baselinestartind != nothing && baselineendind != nothing)
        baselinestrarray = strip.(chomp.(transstrarray[(baselinestartind+
                                                        1):(baselineendind-1)]));
        if (length(baselinestrarray) > 0)
            for currbaselinestr in baselinestrarray
                currbaselinestrsplit = split(currbaselinestr);
                currmollabel = currbaselinestrsplit[1];
                currmolind = findfirst(myMolData.molLabels .== currmollabel);
                
                if (currmolind != nothing)
                    ## Note: all Euler angles should be given in the transformation file in DEGREES
                    baselinetrans[3*myMolData.nummols+3*currmolind-2] = parse(FT, currbaselinestrsplit[2]);
                    baselinetrans[3*myMolData.nummols+3*currmolind-1] = parse(FT, currbaselinestrsplit[3]);
                    baselinetrans[3*myMolData.nummols+3*currmolind] = parse(FT, currbaselinestrsplit[4]);
                    if (length(currbaselinestrsplit) == 7)
                        baselinetrans[3*currmolind-2] = parse(FT, currbaselinestrsplit[5]);
                        baselinetrans[3*currmolind-1] = parse(FT, currbaselinestrsplit[6]);
                        baselinetrans[3*currmolind] = parse(FT, currbaselinestrsplit[7]);
                    end
                end
            end
        end
    end

    # sort such that the displacement of the last molecule changes first
    # between transformation labels, then the displacement of the
    # next-to-last molecule, and so on, then the rotation of the last
    # molecule, and so on
    allmoltranssorted = sortslices(allmoltrans, dims=2);
    sortindvec = round.(Int, vec(allmoltranssorted[end, :]));
    permute!(allmoltranslabels, sortindvec);
    allmoltrans = allmoltranssorted[1:end-1, :];

    # convert Euler angles into (alibi) rotation matrices
    allmoltranslations = deepcopy(allmoltrans[3*myMolData.nummols+1:end, :]);
    allmolrotations = Array{FT, 2}(undef, 3*myMolData.nummols, 3*numtrans);
    for tt=1:numtrans
        for nn=1:myMolData.nummols
            if (allmoltrans[3*nn-2:3*nn, tt] == zeros(FT, 3))
                allmolrotations[3*nn-2:3*nn, 3*tt-2:3*tt] = Matrix{FT}(I, 3, 3);
            else
                cosphi = cosd(allmoltrans[3*nn-2, tt]);
                sinphi = sind(allmoltrans[3*nn-2, tt]);
                costheta = cosd(allmoltrans[3*nn-1, tt]);
                sintheta = sind(allmoltrans[3*nn-1, tt]);
                cospsi = cosd(allmoltrans[3*nn, tt]);
                sinpsi = sind(allmoltrans[3*nn, tt]);
                
                allmolrotations[3*nn-2, 3*tt-2] = costheta*cospsi;
                allmolrotations[3*nn-2, 3*tt-1] = costheta*sinpsi;
                allmolrotations[3*nn-2, 3*tt] = -1.0*sintheta;
                allmolrotations[3*nn-1, 3*tt-2] = sinphi*sintheta*cospsi - cosphi*sinpsi;
                allmolrotations[3*nn-1, 3*tt-1] = sinphi*sintheta*sinpsi + cosphi*cospsi;
                allmolrotations[3*nn-1, 3*tt] = sinphi*costheta;
                allmolrotations[3*nn, 3*tt-2] = cosphi*sintheta*cospsi + sinphi*sinpsi;
                allmolrotations[3*nn, 3*tt-1] = cosphi*sintheta*sinpsi - sinphi*cospsi;
                allmolrotations[3*nn, 3*tt] = cosphi*costheta;
            end
        end
    end
    
    baselinetranslations = deepcopy(baselinetrans[3*myMolData.nummols+1:end]);
    baselinerotations = Array{FT, 2}(undef, 3*myMolData.nummols, 3);
    for nn=1:myMolData.nummols
        if (baselinetrans[3*nn-2:3*nn] == zeros(FT, 3))
            baselinerotations[3*nn-2:3*nn, :] = Matrix{FT}(I, 3, 3);
        else
            cosphi = cosd(baselinetrans[3*nn-2]);
            sinphi = sind(baselinetrans[3*nn-2]);
            costheta = cosd(baselinetrans[3*nn-1]);
            sintheta = sind(baselinetrans[3*nn-1]);
            cospsi = cosd(baselinetrans[3*nn]);
            sinpsi = sind(baselinetrans[3*nn]);
            
            baselinerotations[3*nn-2, 1] = costheta*cospsi;
            baselinerotations[3*nn-2, 2] = costheta*sinpsi;
            baselinerotations[3*nn-2, 3] = -1.0*sintheta;
            baselinerotations[3*nn-1, 1] = sinphi*sintheta*cospsi - cosphi*sinpsi;
            baselinerotations[3*nn-1, 2] = sinphi*sintheta*sinpsi + cosphi*cospsi;
            baselinerotations[3*nn-1, 3] = sinphi*costheta;
            baselinerotations[3*nn, 1] = cosphi*sintheta*cospsi + sinphi*sinpsi;
            baselinerotations[3*nn, 2] = cosphi*sintheta*sinpsi - sinphi*cospsi;
            baselinerotations[3*nn, 3] = cosphi*costheta;
        end
    end

    # after sorting, determine if each successive transformation affects a given molecule
    changedfrombefore = trues(myMolData.nummols, numtrans);
    for tt=2:numtrans
        for currmolind=1:myMolData.nummols
            if (allmoltrans[3*myMolData.nummols+3*currmolind-2, tt] ==
                allmoltrans[3*myMolData.nummols+3*currmolind-2, tt-1] &&
                allmoltrans[3*myMolData.nummols+3*currmolind-1, tt] ==
                allmoltrans[3*myMolData.nummols+3*currmolind-1, tt-1] &&
                allmoltrans[3*myMolData.nummols+3*currmolind, tt] ==
                allmoltrans[3*myMolData.nummols+3*currmolind, tt-1] &&
                allmoltrans[3*currmolind-2, tt] == allmoltrans[3*currmolind-2, tt-1] &&
                allmoltrans[3*currmolind-1, tt] == allmoltrans[3*currmolind-1, tt-1] &&
                allmoltrans[3*currmolind, tt] == allmoltrans[3*currmolind, tt-1])

                changedfrombefore[currmolind, tt] = false;

            end
        end
    end

    # horizontally concatenate all atomic coordinates over all (not just unique) molecules
    allAtomPos = Array{FT, 2}(undef, 3, myMolData.cumnumatomslist[end]);

    allAtomPos[:, 1:myMolData.numatomslist[1]] .=
    myMolData.uniqueMolDataArray[1].atompos0[:, :];
    
    for nn=2:myMolData.nummols
      
        currmolDupeInd = myMolData.dupemolinds[nn];
        if (currmolDupeInd == -1)
            allAtomPos[:,
                       myMolData.cumnumatomslist[nn-1]+1:myMolData.cumnumatomslist[nn]] .=
            myMolData.uniqueMolDataArray[myMolData.uniquemolinds[nn]].atompos0[:, :];
        else
            allAtomPos[:,
                       myMolData.cumnumatomslist[nn-1]+1:myMolData.cumnumatomslist[nn]] .=
            myMolData.uniqueMolDataArray[myMolData.dupemolinds[nn]].atompos0[:, :];
        end

    end

    # create and return the TransData
    myTransData = TransData{FT}(allmoltranslations, allmolrotations,
                                allmoltranslabels, changedfrombefore,
                                numtrans, allAtomPos, baselinetranslations,
                                baselinerotations);
    initializeAtomPos!(myTransData, myMolData);
    return myTransData;
    
end

"""

    initializeAtomPos!(myTransData::TransData{FT},
                       myMolData::MolSystem{<:OneMol{FT}}) where FT<:AbstractFloat


Initialize atomic coordinates across all molecules by first shifting coordinates of all
molecules to have centers of mass of 0, then rotating all molecular coordinates according to
baseline orientations of each, then translating to desired baseline centers of mass,
modifying `myTransData` in-place.

See also: [`transformAtomPos!`](@ref), [`untransformAtomPos!`](@ref)

"""
function initializeAtomPos!(myTransData::TransData{FT},
                            myMolData::MolSystem{<:OneMol{FT}}) where FT<:AbstractFloat
    startatom = 1;
    endatom = myMolData.numatomslist[1];
    for nn=1:myMolData.nummols
        if (nn > 1)
            startatom = myMolData.cumnumatomslist[nn-1] + 1;
            endatom = myMolData.cumnumatomslist[nn];
        end
        makeCOM0!(myTransData.allAtomPos, startatom, endatom);
        if (myTransData.baselineRot[3*nn-2:3*nn, :] != Matrix{FT}(I, 3, 3))
            rotatecoords!(myTransData.allAtomPos, myTransData.baselineRot[3*nn-2:3*nn, :],
                          startatom, endatom);
        end
        if (myTransData.baselineTrans[3*nn-2:3*nn] != zeros(FT, 3))
            translateCOM!(myTransData.allAtomPos, myTransData.baselineTrans[3*nn-2:3*nn],
                          startatom, endatom);
        end
    end
end

"""

    transformAtomPos!(myTransData::TransData{FT},
                      myMolData::MolSystem{<:OneMol{FT}},
                      tt::Integer) where FT<:AbstractFloat


Transform atomic coordinates across all molecules by first shifting coordinates of all
molecules to have centers of mass of 0, then rotating all molecular coordinates according to
orientations corresponding to transformation `tt` relative to the baseline, then translating
to desired centers of mass corresponding to transformation `tt` relative to the baseline,
modifying `myTransData` in-place.

See also: [`initializeAtomPos!`](@ref), [`untransformAtomPos!`](@ref)

"""
function transformAtomPos!(myTransData::TransData{FT},
                           myMolData::MolSystem{<:OneMol{FT}},
                           tt::Integer) where FT<:AbstractFloat
    startatom = 1;
    endatom = myMolData.numatomslist[1];
    for nn=1:myMolData.nummols
        if (nn > 1)
            startatom = myMolData.cumnumatomslist[nn-1] + 1;
            endatom = myMolData.cumnumatomslist[nn];
        end
        if (myTransData.baselineTrans[3*nn-2:3*nn] != zeros(FT, 3))
            makeCOM0!(myTransData.allAtomPos, startatom, endatom);
        end
        if (myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] != Matrix{FT}(I, 3, 3))
            rotatecoords!(myTransData.allAtomPos,
                          myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt],
                          startatom, endatom);
        end
        if (myTransData.baselineTrans[3*nn-2:3*nn] .+
            myTransData.transArray[3*nn-2:3*nn, tt] != zeros(FT, 3))
            translateCOM!(myTransData.allAtomPos,
                          myTransData.baselineTrans[3*nn-2:3*nn] .+
                          myTransData.transArray[3*nn-2:3*nn, tt],
                          startatom, endatom);
        end
    end
end
  
"""

    untransformAtomPos!(myTransData::TransData{FT},
                        myMolData::MolSystem{<:OneMol{FT}},
                        tt::Integer) where FT<:AbstractFloat


Undo transformation of atomic coordinates across all molecules by first shifting coordinates
of all molecules to have centers of mass of 0, then rotating all molecular coordinates
according to orientations corresponding to the opposites of transformation `tt` to return to
the baseline, then translating back to the desired baseline centers of mass, modifying
`myTransData` in-place.

See also: [`initializeAtomPos!`](@ref), [`transformAtomPos!`](@ref)

"""
function untransformAtomPos!(myTransData::TransData{FT},
                             myMolData::MolSystem{<:OneMol{FT}},
                             tt::Integer) where FT<:AbstractFloat
    startatom = 1;
    endatom = myMolData.numatomslist[1];
    for nn=1:myMolData.nummols
        if (nn > 1)
            startatom = myMolData.cumnumatomslist[nn-1] + 1;
            endatom = myMolData.cumnumatomslist[nn];
        end
        if (myTransData.baselineTrans[3*nn-2:3*nn] .+
            myTransData.transArray[3*nn-2:3*nn, tt] != zeros(FT, 3))
            makeCOM0!(myTransData.allAtomPos, startatom, endatom);
        end
        if (myTransData.rotArray[3*nn-2:3*nn, 3*tt-2:3*tt] != Matrix{FT}(I, 3, 3))
            rotatecoords!(myTransData.allAtomPos,
                          transpose(myTransData.rotArray[3*nn-2:3*nn,3*tt-2:3*tt]),
                          startatom, endatom);
        end
        if (myTransData.baselineTrans[3*nn-2:3*nn] != zeros(FT, 3))
            translateCOM!(myTransData.allAtomPos,
                          myTransData.baselineTrans[3*nn-2:3*nn],
                          startatom, endatom);
        end
    end
end
