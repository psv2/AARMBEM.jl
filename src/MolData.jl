"""

    OneMol{FT<:AbstractFloat}

Represent properties of a single molecule, with or without phonons, as an abstract type.

See also: [`OneMol_WithPhonons{FT}`](@ref), [`OneMol_NoPhonons{FT}`](@ref)

"""
abstract type OneMol{FT<:AbstractFloat} end;


"""

    OneMol_WithPhonons{FT} <: OneMol{FT}

Represent the properties of a single molecule of ``N`` atoms, with phonons.

`numatoms` is an integer equal to the number of atoms ``N``.

`Qe`, `Me`, `Be`, `Ke`, `MI`, and `BI` are all `numatoms`-element real vectors.

`KI` is a `3*numatoms`-by-`3*numatoms` real matrix.

`molpos0` is a 3-by-`numatoms` matrix where `molpos0[:, p]` represents the coordinates of
the atom labeled ``p``.

`molCOM` is a 3-element real vector containing the unweighted molecular center of mass.

`alphainv` and `GvacGGinf` are `3*numatoms`-by-`3*numatoms` complex matrices respectively
representing the inverse susceptibility and Green's function interaction matrix of the
molecule in isolation in vacuum (allocated but uninitialized until a frequency and
wavevector are specified).

`alpha0` is a `numatoms`-element real vector of the atomic polarizabilities (allocated but
uninitialized until a frequency and wavevector are specified).

`alphae` is an `numatoms`-element complex vector of the purely electronic contributions to
the polarizabilities, useful mainly for vdW computations (allocated but uninitialized until
a frequency and wavevector are specified).

`numblocks1` & `numblocks2` are the integer numbers of neighboring unit cells that have
nontrivial internuclear couplings to the central unit cell (each is 1 molecules that are not
periodic in that direction or have no internuclear couplings/covalent bonds with neighboring
unit cells in that direction).

`blocklist1` & `blocklist2` are range vectors indexing each unit cell block (each has the
single element 0 for molecules that are not periodic in that direction or have no
internuclear couplings/covalent bonds with neighboring unit cells in that direction).

See also: [`OneMol{FT}`](@ref), [`OneMol_NoPhonons{FT}`](@ref) 

"""
struct OneMol_WithPhonons{FT} <: OneMol{FT}

    Qe::Array{FT, 1}
    Me::Array{FT, 1}
    Be::Array{FT, 1}
    Ke::Array{FT, 1}
    MI::Array{FT, 1}
    BI::Array{FT, 1}
    KI::Array{FT, 2}
    
    atompos0::Array{FT, 2}
    numatoms::Integer
    
    alphainv::Array{Complex{FT}, 2}
    GvacGGinf::SharedArray{Complex{FT}, 2}
    GvacGGinfdiag::Array{Complex{FT}, 1}
    alpha0::Array{FT, 1}
    alphae::Array{Complex{FT}, 1}
    
    numblocks1::Integer
    blocklist1::UnitRange{<:Integer}
    numblocks2::Integer
    blocklist2::UnitRange{<:Integer}

end

"""

    OneMol_NoPhonons{FT} <: OneMol{FT}

Represent the properties of a single molecule of ``N`` atoms, without phonons.

`numatoms` is an integer equal to the number of atoms ``N``.

`Qe`, `Me`, `Be`, and `Ke` are all `numatoms`-element real vectors.

`molpos0` is a 3-by-`numatoms` matrix where `molpos0[:, p]` represents the coordinates of
atom ``p``.

`molCOM` is a 3-element real vector containing the unweighted molecular center of mass.

`alphainv` is a `numatoms`-element complex vector representing the local inverse
susceptibility of each atom (allocated but uninitialized until a frequency and wavevector
are specified).

`GvacGGinf` is a `3*numatoms`-by-`3*numatoms` complex matrix representing the Green's
function interaction matrix of the molecule in isolation in vacuum (allocated but
uninitialized until a frequency and wavevector are specified).

`alphae` is `1 ./ alphainv`, useful mainly for vdW computations (allocated but uninitialized
until a frequency and wavevector are specified).

`alpha0` is `abs.(alphae)`, representing real atomic polarizabilities (allocated but
uninitialized until a frequency and wavevector are specified).

See also: [`OneMol{FT}`](@ref), [`OneMol_WithPhonons{FT}`](@ref) 

"""
struct OneMol_NoPhonons{FT} <: OneMol{FT}

    Qe::Array{FT, 1}
    Me::Array{FT, 1}
    Be::Array{FT, 1}
    Ke::Array{FT, 1}
    
    atompos0::Array{FT, 2}
    numatoms::Integer
    
    alphainv::Array{Complex{FT}, 1}
    GvacGGinf::SharedArray{Complex{FT}, 2}
    GvacGGinfdiag::Array{Complex{FT}, 1}
    alpha0::Array{FT, 1}
    alphae::Array{Complex{FT}, 1}

end

"""

    MolSystem{T<:OneMol{<:AbstractFloat}}

Represent the properties of a collection of molecules.

`uniqueMolDataArray` contains data for each unique molecule.

`nummols` is the total number of molecules (unique or not).

`numatomslist` is the list of atoms in every molecule (unique or not).

`cumnumatomslist` is the cumulative sum of the number of atoms in every molecule, to keep
track of overall indexing.

`dupemolinds` is the list of indices tracking duplicate molecules: -1 indicates it is the
first instance of that molecule, while positive index `nn` (in the range `1:nummols`) points
to `nn` as the molecule to be duplicated.

`uniquemolinds` is complementary to `dupemolinds`, listing the index of each unique molecule
in the set of unique molecules, or -1 if it is a duplicate molecule.

See also: [`OneMol{FT}`](@ref), [`OneMol_WithPhonons{FT}`](@ref),
[`OneMol_NoPhonons{FT}`](@ref) 

"""
struct MolSystem{T<:OneMol{<:AbstractFloat}}

    uniqueMolDataArray::Array{T, 1}
    molLabels::Array{<:AbstractString, 1}
  
    nummols::Integer
    numatomslist::Array{<:Integer, 1}
    cumnumatomslist::Array{<:Integer, 1}
    dupemolinds::Array{<:Integer, 1}
    uniquemolinds::Array{<:Integer, 1}

end

function extractKIconst!(KIconst::Array{<:AbstractFloat, 2},
                         molstrarray::Array{<:AbstractString, 1},
                         mollistfilename::AbstractString)

    if (size(KIconst, 1) != 3 || size(KIconst, 2) != 3)
        error("extractKIconst!(): KIconst must be a 3-by-3 matrix");
    end
    
    fill!(KIconst, zero(eltype(KIconst)));
    curridx = findfirst(map(u -> startswith(u, "KIconst "), molstrarray));
    if (curridx != nothing)
        currKIconst = map(u -> parse(eltype(KIconst), u),
                          split(replace(molstrarray[curridx], "KIconst " => "")));
        currKIconstlen = length(currKIconst);
        if (currKIconstlen == 1)
            for ii=1:3
                KIconst[ii, ii] = abs(currKIconst[1]);
            end
        elseif (currKIconstlen == 3)
            for ii=1:3
                KIconst[ii, ii] = abs(currKIconst[ii]);
            end
        elseif (currKIconstlen == 6 || currKIconstlen == 9)
            tempKIconst = zero(KIconst);
            if (currKIconstlen == 6)
        
                tempKIconst[1, 1] = currKIconst[1];
                tempKIconst[2, 1] = currKIconst[2];
                tempKIconst[3, 1] = currKIconst[3];
                tempKIconst[1, 2] = currKIconst[2];
                tempKIconst[2, 2] = currKIconst[4];
                tempKIconst[3, 2] = currKIconst[5];
                tempKIconst[1, 3] = currKIconst[3];
                tempKIconst[2, 3] = currKIconst[5];
                tempKIconst[3, 3] = currKIconst[6];
                
            else
                for ii=1:9
                    tempKIconst[ii] = currKIconst[ii];
                end
            end
            curreigvals, curreigvecs = eigen(Symmetric(tempKIconst));
            KIconst .= curreigvecs * broadcast(*, transpose(curreigvecs),
                                               abs(curreigvals));
        else
            error(string(mollistfilename, ": expected scalar (1)",
                         " or tensor (3, 6, or 9) elements for",
                         " KIconst, got something else"));
        end
    else
        println(string(mollistfilename, ": should optionally",
                       " contain a line of the form: KIconst",
                       " \${KIconstvals...}"));
    end

end

function extractRequiredMolDataVec(datalabel::AbstractString,
                                   molstrarray::Array{<:AbstractString, 1},
                                   mollistfilename::AbstractString, T::DataType=Float64)

    curridx = findfirst(map(u -> startswith(u, string(datalabel, " ")), molstrarray));
    if (curridx == nothing)
        error(string(mollistfilename, ": must contain a line",
                     " of the form: ", datalabel, " \${", datalabel, "filename}"));
    end
    return vec(collect(readdlm(replace(molstrarray[curridx],
               string(datalabel, " ") => ""), ',', T)));

end

function extractOptionalMolDataVec(datalabel::AbstractString, numatoms::Integer,
                                   molstrarray::Array{<:AbstractString, 1},
                                   mollistfilename::AbstractString, T::DataType=Float64)

    datavec = zeros(T, numatoms);
    curridx = findfirst(map(u -> startswith(u, string(datalabel, " ")), molstrarray));
    if (curridx != nothing)
        datavec = vec(collect(readdlm(replace(molstrarray[curridx],
                      string(datalabel, " ") => ""), ',', T)))
    else
        println(string(mollistfilename, ": should optionally, as",
                       " recommended, contain a line of the form: ", datalabel,
                       " \${", datalabel, "filename}"));
    end
    return datavec;

end

function extractRequiredMolDataMat(datalabel::AbstractString,
                                   molstrarray::Array{<:AbstractString, 1},
                                   mollistfilename::AbstractString, T::DataType=Float64)

    curridx = findfirst(map(u -> startswith(u, string(datalabel, " ")), molstrarray));
    if (curridx == nothing)
        error(string(mollistfilename, ": must contain a line",
                     " of the form: ", datalabel, " \${", datalabel, "filename}"));
    end
    return readdlm(replace(molstrarray[curridx],
                   string(datalabel, " ") => ""), ',', T);

end

function extractBescale(molstrarray::Array{<:AbstractString, 1},
                        mollistfilename::AbstractString, T::DataType=Float64)

    curridx = findfirst(map(u -> startswith(u, "Bescale "), molstrarray));
    Bescale = one(T);
    if (curridx != nothing)
        Bescale = parse(T, replace(molstrarray[curridx], "Bescale " => ""));
    else
        println(string(mollistfilename, ": should optionally",
                       " contain a line of the form: Bescale \${Bescaleval}"));
    end
    return Bescale;

end

function extractKIconstatomchoices(molstrarray::Array{<:AbstractString, 1},
                                   mollistfilename::AbstractString, numatoms::Integer)

    curridx1 = findfirst(map(u -> startswith(u, "KIconstatom1 "), molstrarray));
    curridx2 = findfirst(map(u -> startswith(u, "KIconstatom2 "), molstrarray));
    atom1 = 1;
    atom2 = numatoms;

    if (curridx1 != nothing)
        atom1 = parse(Int, replace(molstrarray[curridx1], "KIconstatom1 " => ""));
    else
        println(string(mollistfilename, ": should optionally",
                       " contain a line of the form: KIconstatom1 \${KIconstatom1}"));
        println("(default: KIconstatom1 is set to 1)");
    end
    if (curridx2 != nothing)
        atom2 = parse(Int, replace(molstrarray[curridx2], "KIconstatom2 " => ""));
    else
        println(string(mollistfilename, ": should optionally",
                       " contain a line of the form: KIconstatom2 \${KIconstatom2}"));
        println("(default: KIconstatom2 is set to numatoms)");
    end

    if (atom1 < 1 || atom1 > numatoms)
        atom1 = 1;
        println(string(mollistfilename, ": required range is ",
                       "1 <= KIconstatom1 < KIconstatom2 <= numatoms"));
        println("(KIconstatom1 has been reset to 1)");
    end
    if (atom2 < 1 || atom2 > numatoms)
        atom2 = numatoms;
        println(string(mollistfilename, ": required range is ",
                       "1 <= KIconstatom1 < KIconstatom2 <= numatoms"));
        println("(KIconstatom2 has been reset to numatoms)");
    end
    if (atom2 < atom1)
        tempatom = atom2;
        atom2 = atom1;
        atom1 = tempatom;
        println(string(mollistfilename, ": recommended order is ",
                       "1 <= KIconstatom1 < KIconstatom2 <= numatoms"));
        println("(KIconstatom1 and KIconstatom2 have been switched)");
    end
    if (atom1 == atom2)
      println(string(mollistfilename, ": required order being ",
                     "1 <= KIconstatom1 < KIconstatom2 <= numatoms ",
                     "means KIconstatom1 and KIconstatom2 should be distinct"));
        if (atom2 == 1)
            atom2 = numatoms;
            println("(KIconstatom2 has been reset to numatoms)");
        else
            atom1 = 1;
            println("(KIconstatom1 has been reset to 1)");
        end
    end

    return atom1, atom2;

end

function extractGeneralMolData(molstrarray::Array{<:AbstractString, 1},
                               mollistfilename::AbstractString, T::DataType=Float64)

    Qe = extractRequiredMolDataVec("Qe", molstrarray, mollistfilename, T);
    Me = extractRequiredMolDataVec("Me", molstrarray, mollistfilename, T);
    Ke = extractRequiredMolDataVec("Ke", molstrarray, mollistfilename, T);
    numatoms = length(Qe);
    Be = extractOptionalMolDataVec("Be", numatoms, molstrarray, mollistfilename, T);
    Bescale = extractBescale(molstrarray, mollistfilename, T);
    lmul!(Bescale, Be);
    xyz = extractRequiredMolDataMat("xyz", molstrarray, mollistfilename, T);
    return numatoms, Qe, Me, Be, Ke, xyz;

end

function addKIconsttoKIselectedatoms!(KI::Array{FT, 2}, KIconst::Array{FT, 2},
                                      numatoms::Integer, atom1::Integer,
                                      atom2::Integer) where FT<:AbstractFloat

    if (size(KIconst, 1) != 3 || size(KIconst, 2) != 3)
        error("KIconst must be a 3-by-3 array");
    end
    if (rem(size(KI, 1), 3*numatoms) != 0 || rem(size(KI, 2), 3*numatoms) != 0)
        error("KI must have each dimension be an integer multiple of 3N,",
              " where N is numatoms");
    end

    numblocks1 = div(size(KI, 1), 3*numatoms);
    numblocks2 = div(size(KI, 2), 3*numatoms);

    if (iseven(numblocks1) || iseven(numblocks2))
        error("KI must have each dimension be an integer multiple of 3N,",
              " where N is numatoms, and where each multiple must be odd");
    end
    
    blocklist1 = div(1 - numblocks1, 2):div(numblocks1 - 1, 2);
    blocklist2 = div(1 - numblocks2, 2):div(numblocks2 - 1, 2);
    if (numblocks1 == 1 && numblocks2 == 1)
        KI[3*atom1-2:3*atom1, 3*atom1-2:3*atom1] .+= KIconst[:, :];
        KI[3*atom2-2:3*atom2, 3*atom2-2:3*atom2] .+= KIconst[:, :];
    else
        KI[3*numatoms*div(numblocks1 - 1, 2)+3*atom1-2:3*numatoms*div(numblocks1 - 1,
                                                                      2)+3*atom1,
           3*numatoms*div(numblocks2 - 1, 2)+3*atom1-2:3*numatoms*div(numblocks2 - 1,
                                                                      2)+3*atom1] .+=
           KIconst[:, :];
    end
    
    return numblocks1, numblocks2, blocklist1, blocklist2;

end

"""

    read1MolDataWithPhonons(molstrarray::Array{<:AbstractString, 1},
                            mollistfilename::AbstractString, T::DataType=Float64)

Read a section of a configuration file `molstrarray` (an array of strings) corresponding to
molecular data and return the data corresponding to that molecule, with phonons, of type
`OneMol_WithPhonons{T}`. Create data types to use `T` (satisfying `T<:AbstractFloat`) for
storing floating-point numbers, and output any errors with reference to the configuration file
`mollistfilename` (a string).

See also: [`OneMol{FT}`](@ref), [`OneMol_WithPhonons{FT}`](@ref) 

"""
function read1MolDataWithPhonons(molstrarray::Array{<:AbstractString, 1},
                                 mollistfilename::AbstractString, T::DataType=Float64)

    if (!(T <: AbstractFloat) || isabstracttype(T))
        error("read1MolDataWithPhonons(): T must be a concrete subtype of AbstractFloat");
    end
  
    KIconst = zeros(T, 3, 3);
    extractKIconst!(KIconst, molstrarray, mollistfilename);
    
    numatoms, Qe, Me, Be, Ke, molpos0 = extractGeneralMolData(molstrarray,
                                                              mollistfilename, T);
    
    MI = extractRequiredMolDataVec("MI", molstrarray, mollistfilename, T);
    BI = extractOptionalMolDataVec("BI", numatoms, molstrarray, mollistfilename, T);

    atom1, atom2 = extractKIconstatomchoices(molstrarray, mollistfilename, numatoms);
    
    KI = extractRequiredMolDataMat("KI", molstrarray, mollistfilename, T);
    
    numblocks1, numblocks2, blocklist1, blocklist2 =
      addKIconsttoKIselectedatoms!(KI, KIconst, numatoms, atom1, atom2);

    GvacGG = SharedArray{Complex{T}, 2}(3*numatoms, 3*numatoms);
    fill!(GvacGG, zero(eltype(GvacGG)));
    
    return OneMol_WithPhonons{T}(Qe, Me, Be, Ke, MI, BI, KI, molpos0,
                                 numatoms, zeros(Complex{T}, 3*numatoms, 3*numatoms),
                                 GvacGG, zeros(Complex{T}, numatoms), zeros(T, numatoms),
                                 zeros(Complex{T}, numatoms), numblocks1, blocklist1,
                                 numblocks2, blocklist2);
end

"""

    read1MolDataNoPhonons(molstrarray::Array{<:AbstractString, 1},
                          mollistfilename::AbstractString, T::DataType=Float64)

Read a section of a configuration file `molstrarray` (an array of strings) corresponding to
molecular data and return the data corresponding to that molecule, without phonons, of type
`OneMol_NoPhonons{T}`. Create data types to use `T` (satisfying `T<:AbstractFloat`) for
storing floating-point numbers, and output any errors with reference to the configuration
file `mollistfilename` (a string).

See also: [`OneMol{FT}`](@ref), [`OneMol_NoPhonons{FT}`](@ref) 

"""
function read1MolDataNoPhonons(molstrarray::Array{<:AbstractString, 1},
                               mollistfilename::AbstractString, T::DataType=Float64)

    if (!(T<:AbstractFloat) || isabstracttype(T))
        error("read1MolDataWithPhonons(): T must be a concrete subtype of AbstractFloat");
    end

    numatoms, Qe, Me, Be, Ke, molpos0 = extractGeneralMolData(molstrarray,
                                                              mollistfilename, T);
    
    GvacGG = SharedArray{Complex{T}, 2}(3*numatoms, 3*numatoms);
    fill!(GvacGG, zero(eltype(GvacGG)));
    
    return OneMol_NoPhonons{T}(Qe, Me, Be, Ke, molpos0, numatoms,
                               zeros(Complex{T}, numatoms), GvacGG,
                               zeros(Complex{T}, numatoms), zeros(T, numatoms),
                               zeros(Complex{T}, numatoms));
end

"""

    readMolSystem(mollistfilename::AbstractString,
                  nophonons:Bool, T::DataType=Float64)

Read a configuration file `mollistfilename` (a string) to extract molecular data of type
either `OneMol_NoPhonons{T}`, if `nophonons` is `true`, or `OneMol_WithPhonons{T}`
otherwise, returning information about the composite system of type `MolSystem{molType{T}}`.
Ensure data types use `T` (satisfying `T<:AbstractFloat`) for storing floating-point
numbers.

See also: [`OneMol{T}`](@ref), [`OneMol_WithPhonons{T}`](@ref),
[`OneMol_NoPhonons{T}`](@ref), [`MolSystem{molType{T}}`](@ref)

"""
function readMolSystem(mollistfilename::AbstractString,
                       nophonons::Bool, T::DataType=Float64)

    if (!(T<:AbstractFloat) || isabstracttype(T))
        error("read1MolDataWithPhonons(): T must be a concrete subtype of AbstractFloat");
    end
  
    myreadfunc = read1MolDataNoPhonons;
    molType = OneMol_NoPhonons;
    if (!nophonons)
        molType = OneMol_WithPhonons;
        myreadfunc = read1MolDataWithPhonons;
    end
  
    molstrarray = strip.(chomp.(readlines(mollistfilename)));
    molstartinds = findall(map(u -> startswith(u, "MOL "), molstrarray));
    molendinds = findall(molstrarray .== "ENDMOL");
    nummols = length(molendinds);
    if (length(molstartinds) != nummols)
        error(string(mollistfilename, ": number of MOL and ENDMOL statements",
                     " should be the same"));
    end
    numatomslist = zeros(Int, nummols);
    molbothinds = vec(transpose(hcat(molstartinds, molendinds)));
    if (molbothinds != sort(molbothinds))
        error(string(mollistfilename, ": MOL and ENDMOL statements",
                     " should alternate in ordered occurrence"));
    end

    molstrarraydupeinds = findall(map(u -> startswith(u, "DUPLICATE "), molstrarray));
    numdupes = length(molstrarraydupeinds);
    numunique = nummols - numdupes;
    dupemolinds = fill(-1, nummols);
    
    uniqueMolDataArray = Array{molType{T}, 1}(undef, numunique);
    allMolLabelArray = Array{String, 1}(undef, nummols);
    uniquemolcounter = zero(Int);
    uniquemolinds = fill(-1, nummols);
    for nn=1:nummols
        allMolLabelArray[nn] = replace(strip(chomp(molstrarray[molstartinds[nn]])),
                                       "MOL " => "");
        currmolstrarray = strip.(chomp.(molstrarray[molstartinds[nn]+1:molendinds[nn]-1]));
        if (findfirst(map(u -> startswith(u, "DUPLICATE "), currmolstrarray)) == nothing)
            uniquemolcounter += 1;
            uniqueMolDataArray[uniquemolcounter] = myreadfunc(currmolstrarray,
                                                              mollistfilename, T);
            numatomslist[nn] = uniqueMolDataArray[uniquemolcounter].numatoms;
            uniquemolinds[nn] = uniquemolcounter;
        elseif (nn > 1)
            dupemolinds[nn] = findfirst(allMolLabelArray .==
                                        replace(strip(chomp(currmolstrarray[1])),
                                                "DUPLICATE " => ""));
            numatomslist[nn] = uniqueMolDataArray[dupemolinds[nn]].numatoms;
        else
            error(string(mollistfilename, ": the first molecule listed cannot be a",
                         "duplicate"));
        end
    end
    cumnumatomslist = cumsum(numatomslist);

    return MolSystem{molType{T}}(uniqueMolDataArray, allMolLabelArray,
                                 nummols, numatomslist, cumnumatomslist,
                                 dupemolinds, uniquemolinds);

end
