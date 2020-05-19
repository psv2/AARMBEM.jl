"""

    readRequiredFilenameARG(myARGS, currfilename, callingfunc)

Read from the list of command-line arguments `myARGS` the file name prefixed by
`currfilename` and an equals sign. Throw an error (using `callingfunc` to emphasize where it
is needed) if that argument is missing.

See also: [`readOptionalFilenameARG`](@ref), [`readGenvstr`](@ref),
[`readOptionalBoolcondARG`](@ref), [`readOptionalFloatType`](@ref),
[`readfreqlist!`](@ref), [`readklist`](@ref)

"""
function readRequiredFilenameARG(myARGS, currfilename::AbstractString,
                                 callingfunc::AbstractString)

    curridx = findfirst(map(u -> startswith(u, string(currfilename, "=")), myARGS));
    if (curridx == nothing)
        error(string("In", callingfunc, ", ", currfilename, "=\${", currfilename, "} is a ",
                     " REQUIRED argument"));
    end
    return replace(myARGS[curridx], string(currfilename, "=") => "");

end

"""

    readOptionalFilenameARG(myARGS, currfilename)

Read from the list of command-line arguments `myARGS` the file name prefixed by
`currfilename` and an equals sign. Return a blank string if that argument is missing.

See also: [`readRequiredFilenameARG`](@ref), [`readGenvstr`](@ref),
[`readOptionalBoolcondARG`](@ref), [`readOptionalFloatType`](@ref),
[`readfreqlist!`](@ref), [`readklist`](@ref)

"""
function readOptionalFilenameARG(myARGS, currfilename::AbstractString)

    curridx = findfirst(map(u -> startswith(u, string(currfilename, "=")), myARGS));
    if (curridx != nothing)
        return replace(myARGS[curridx], string(currfilename, "=") => "");
    end
    return "";

end

"""

    readGenvstr(myARGS)

Read from the list of command-line arguments `myARGS` the environmental Green's function
prefixed by `\"Genv=\"`. Return `\"PEC\"` (representing a PEC plane) only if that is
explicitly found, otherwise return `\"VAC\"` by default.

See also: [`readRequiredFilenameARG`](@ref), [`readOptionalFilenameARG`](@ref),
[`readOptionalBoolcondARG`](@ref), [`readOptionalFloatType`](@ref),
[`readfreqlist!`](@ref), [`readklist`](@ref)

"""
function readGenvstr(myARGS) 
    
    if (findfirst(myARGS .== "Genv=PEC") != nothing) 
        return "PEC";
    end
    return "VAC";
    
end

"""

    readOptionalBoolcondARG(myARGS, currvarname, outfilename)

Read from the list of command-line arguments `myARGS` optional boolean arguments labeled by
`currvarname`, and modify `outfilename` as appropriate.

See also: [`readRequiredFilenameARG`](@ref), [`readOptionalFilenameARG`](@ref),
[`readGenvstr`](@ref), [`readOptionalFloatType`](@ref),
[`readfreqlist!`](@ref), [`readklist`](@ref)

"""
function readOptionalBoolcondARG(myARGS, currvarname::AbstractString,
                                 outfilename::AbstractString)

    currboolvar = (findfirst(myARGS .== currvarname) != nothing);
    outfilename = currboolvar ? string(currvarname, "_", outfilename) : outfilename;
    return currboolvar, outfilename;

end

"""

    readOptionalFloatType(myARGS)

Read from the list of command-line arguments `myARGS` optional arguments labeling the
floating-point type. Return `Float64` by default if nothing else is found.

See also: [`readRequiredFilenameARG`](@ref), [`readOptionalFilenameARG`](@ref),
[`readGenvstr`](@ref), [`readOptionalBoolcondARG`](@ref),
[`readfreqlist!`](@ref), [`readklist`](@ref)

"""
function readOptionalFloatType(myARGS)

    FT = Float64;
    curridx = findfirst(map(u -> startswith(u, "FloatType="), myARGS));
    if (curridx != nothing)
        currsymb = Symbol(replace(myARGS[curridx], "FloatType=" => ""));
        if (isdefined(Main, currsymb))
            currFT = getfield(Main, currsymb);
            if (isa(currFT, DataType) && currFT<:AbstractFloat)
                FT = currFT;
            end
        end
    end
    return FT;

end

"""

    readfreqlist!(freqlist, myARGS, callingfunc)

Read from the list of command-line arguments `myARGS` the list of frequencies, modifying
`freqlist` \"in-place\" (actually appending as needed), and using the name `callingfunc`
to be more informative when throwing errors.

!!! note

    Any frequencies will be converted to their absolute values. Frequencies are not allowed
    to be exactly zero.

See also: [`readRequiredFilenameARG`](@ref), [`readOptionalFilenameARG`](@ref),
[`readGenvstr`](@ref), [`readOptionalBoolcondARG`](@ref),
[`readOptionalFloatType`](@ref), [`readklist`](@ref)

"""
function readfreqlist!(freqlist::Array{FT, 1}, myARGS,
                       callingfunc::AbstractString) where FT<:AbstractFloat
  
    freq = zero(eltype(freqlist));
    if (findfirst(map(u -> startswith(u, "freqlistfilename="), myARGS)) != nothing)
      
        freqliststrarray = strip.(chomp.(readlines(
                                             replace(myARGS[findfirst(map(u ->
                                                                  startswith(u, "freqlistfilename="),
                                                                  myARGS))],
                                             "freqlistfilename=" => ""))));
        filter!(x -> !isempty(x), freqliststrarray);
        append!(freqlist, abs.(parse.(Complex{FT}, freqliststrarray)));
        if (findfirst(freqlist .== zero(freqlist)) != nothing)
            error("In ", callingfunc, ", argument freqlistfilename=\${freqlistfilename}",
                  " cannot have any entries that are exactly 0.0");
        end
        print("Warning: In ", callingfunc, ", all entries in argument",
              " freqlistfilename=\${freqlistfilename} are converted to their");
        
    elseif (findfirst(map(u -> startswith(u, "freq="), myARGS)) != nothing)
      
        freq = abs(parse(Complex{FT},
                         replace(myARGS[findfirst(map(u -> startswith(u, "freq="),
                                                      myARGS))], "freq=" => "")));
        if (freq == zero(freq))
            error("In ", callingfunc, ", argument freq=\${freq} cannot be exactly 0.0");
        end
        print("Warning: In ", callingfunc, ", argument freq=\${freq} is converted to its");
        push!(freqlist, freq);
        
    else
        error("In ", callingfunc, ", EITHER arguments freq=\${freq} OR",
              " freqlistfilename=\${freqlistfilename} are REQUIRED");
    end
    print(" absolute value(s)");
    if (callingfunc == "AARMBEM_vdW.jl")
        println(", then multiplied by 1im");
    else
        println("");
    end
        
end

"""

    readklist(myARGS, callingfunc, numdims)

Read from the list of command-line arguments `myARGS` the list of Bloch wavevectors,
returning `zeros(FT, 1, 3)` if `numdims` is zero.

See also: [`readRequiredFilenameARG`](@ref), [`readOptionalFilenameARG`](@ref),
[`readGenvstr`](@ref), [`readOptionalBoolcondARG`](@ref),
[`readOptionalFloatType`](@ref), [`readfreqlist!`](@ref)

"""
function readklist(myARGS, callingfunc::AbstractString, numdims::Integer,
                   FT::DataType=Float64)

    if (!(FT <: AbstractFloat))
        error("FT must be a subtype of AbstractFloat");
    end
    klist = zeros(FT, 1, 3);
    if (numdims <= 0)
      
        return klist;
        
    elseif (findfirst(map(u -> startswith(u, "klistfilename="), myARGS)) != nothing)

        return readdlm(replace(myARGS[findfirst(map(u ->
                                                    startswith(u, "klistfilename="),
                                                    myARGS))],
                               "klistfilename=" => ""), ' ', FT);

    elseif (findfirst(map(u -> startswith(u, "kx="), myARGS)) != nothing ||
                findfirst(map(u -> startswith(u, "ky="), myARGS)) != nothing ||
                findfirst(map(u -> startswith(u, "kz="), myARGS)) != nothing)

        if (findfirst(map(u -> startswith(u, "kx="), myARGS)) != nothing)
                
            klist[1, 1] = parse(FT,
                                replace(myARGS[findfirst(map(u ->
                                                             startswith(u, "kx="),
                                                             myARGS))],
                                        "kx=" => ""));
            
        end
        if (findfirst(map(u -> startswith(u, "ky="), myARGS)) != nothing)
                  
            klist[1, 2] = parse(FT,
                                replace(myARGS[findfirst(map(u ->
                                                             startswith(u, "ky="),
                                                             myARGS))],
                                        "ky=" => ""));
            
        end
        if (findfirst(map(u -> startswith(u, "kz="), myARGS)) != nothing)
          
            klist[1, 3] = parse(FT,
                                replace(myARGS[findfirst(map(u ->
                                                             startswith(u, "kz="),
                                                             myARGS))],
                                        "kz=" => ""));
            
        end

        return klist;
        
    else
        error("In ", callingfunc, ", if there are nontrivial periodic dimensions in",
              " \${periodicfilename}, EITHER arguments",
              " klistfilename=\${klistfilename} OR AT LEAST ONE OF kx=\${kx} OR",
              " ky=\${ky} OR kz=\${kz} are REQUIRED");
    end

end
