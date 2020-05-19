"""

    PeriodicData{FT<:AbstractFloat}

Store information relevant to spatial periodicity in a system. `numdims` should be 0, 1, or
2, representing the number of periodic dimensions. `latticevecs` and `reciprocalvecs` should
each by 3-by-2 matrices, with each column representing a 3-element lattice or corresponding
reciprocal vector, with the vector being zeros along a dimension with no periodicity.

"""
struct PeriodicData{FT<:AbstractFloat}

    numdims::Integer
    latticevecs::Array{FT, 2}
    reciprocalvecs::Array{FT, 2}

end

"""

    readPeriodicData(FT::DataType=Float64, periodicfilename::AbstractString="")

Construct PeriodicData from an appropriate floating-point data type `FT` and a file of name
`periodicfilename` containing information about spatial periodicity; if no file is
specified, produce a default set of periodic data (no spatial periodicity at all).

"""
function readPeriodicData(FT::DataType=Float64, periodicfilename::AbstractString="")

    if (!(FT<:AbstractFloat) || isabstracttype(FT))
        error("readPeriodicData(): FT must be a concrete subtype of AbstractFloat");
    end

    if (periodicfilename == "")
        return PeriodicData(zero(Integer), zeros(FT, 3, 2), zeros(FT, 3, 2));
    end
    latticedata = vec(readdlm(periodicfilename, ',', FT));
    if (length(latticedata) == 3)
        return PeriodicData(one(Integer), hcat(latticedata, zero(latticedata)),
                            hcat(2*pi .* latticedata/dot(latticedata, latticedata),
                                 zero(latticedata)));
    elseif (length(latticedata) == 6)
        latticevecs = reshape(latticedata, 3, 2);
        eperp = cross(latticevecs[:,1], latticevecs[:,2]);
        Auc = norm(eperp);
        eperp = eperp ./ Auc;
        return PeriodicData(2, latticevecs,
                            (2*pi/Auc).*hcat(cross(latticevecs[:,2], eperp),
                                                   cross(eperp, latticevecs[:,1])));
    else
      error(string("readPeriodicData(): ", periodicfilename, " should have exactly 3 (for",
                   " 1D periodicity) or 6 (for 2D periodicity) real elements, each on a",
                   " separate line"));
    end
end

"""

    constructlatticefromreciprocal!(myPeriodicData::PeriodicData{<:AbstractFloat})

Reset the lattice vectors of `myPeriodicData` from its reciprocal vectors.

"""
function constructlatticefromreciprocal!(myPeriodicData::PeriodicData{<:AbstractFloat})

    if (myPeriodicData.numdims == 0)
        fill!(myPeriodicData.latticevecs, zero(eltype(myPeriodicData.latticevecs)));
    elseif (myPeriodicData.numdims == 1)
        myPeriodicData.latticevecs[:, 1] .= 2*pi*myPeriodicData.reciprocalvecs[:, 1] /
        dot(myPeriodicData.reciprocalvecs[:, 1], myPeriodicData.reciprocalvecs[:, 1]);
        myPeriodicData.latticevecs[:, 2] .= zeros(eltype(myPeriodicData.latticevecs), 3);
    elseif (myPeriodicData.numdims == 2)
        eperp = cross(myPeriodicData.reciprocalvecs[:, 1],
                      myPeriodicData.reciprocalvecs[:, 2]);
        ABZ = norm(eperp);
        eperp = eperp ./ ABZ;
        myPeriodicData.latticevecs .= (1/ABZ) .*
        hcat(cross(myPeriodicData.reciprocalvecs[:,2], eperp),
             cross(eperp, myPeriodicData.reciprocalvecs[:,1]));
    else
        error("PeriodicData: numdims must be 0, 1, or 2");
    end
end

"""

    constructreciprocalfromlattice!(myPeriodicData::PeriodicData{<:AbstractFloat})

Reset the reciprocal vectors of `myPeriodicData` from its lattice vectors.

"""
function constructreciprocalfromlattice!(myPeriodicData::PeriodicData{<:AbstractFloat})

    if (myPeriodicData.numdims == 0)
        fill!(myPeriodicData.reciprocalvecs, zero(eltype(myPeriodicData.latticevecs)));
    elseif (myPeriodicData.numdims == 1)
        myPeriodicData.reciprocalvecs[:, 1] .= 2*pi*myPeriodicData.latticevecs[:, 1] /
        dot(myPeriodicData.latticevecs[:, 1], myPeriodicData.latticevecs[:, 1]);
        myPeriodicData.reciprocalvecs[:, 2] .= zeros(eltype(myPeriodicData.reciprocalvecs),
                                                     3);
    elseif (myPeriodicData.numdims == 2)
        eperp = cross(myPeriodicData.latticevecs[:, 1], myPeriodicData.latticevecs[:, 2]);
        Auc = norm(eperp);
        eperp = eperp ./ Auc;
        myPeriodicData.reciprocalvecs .= (2*pi/Auc) .*
        hcat(cross(myPeriodicData.latticevecs[:, 2], eperp),
             cross(eperp, myPeriodicData.latticevecs[:, 1]));
    else
        error("PeriodicData: numdims must be 0, 1, or 2");
    end
end

