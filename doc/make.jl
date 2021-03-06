using Documenter;
push!(LOAD_PATH, "../src/");
using SharedArrays, SpecialFunctions, Distributed, Printf
using DelimitedFiles, LinearAlgebra
include("../src/AARMBEM_ReadCMDARGS.jl");
include("../src/E1func.jl");
include("../src/GenericMatrixVectorUtils.jl");
include("../src/MolData.jl");
include("../src/TransData.jl");
include("../src/PeriodicData.jl");
include("../src/MolTransUtils.jl");
include("../src/MolAlpha.jl");
include("../src/GFVACGG.jl");
include("../src/GFPECGG.jl");
include("../src/MolAssembleGFGG.jl");
include("../src/heattransfer.jl");
include("../src/vdWenergy.jl");

makedocs(
    sitename = "AARMBEM.jl",
    authors = string("Prashanth S. Venkataram, ",
                     "Jan Hermann, ",
                     "Alexandre Tkatchenko, ",
                     "and Alejandro W. Rodriguez"),
    pages = [
        "Home" => "index.md",
        "Background" => [
            "History" => "theory/historicalbackground.md",
            "Theory" => "theory/theoreticalbackground.md",
            "Common acronyms" => "theory/acronyms.md"
        ],
        "Implementation Notes" => [
            "Units" => "implementation/units.md",
            "Matrix storage conventions" => "implementation/matrixstorage.md",
            "Data structures and rigid transformations" => "implementation/datastructuresandtransformations.md",
            "Performance notes and required packages" => "implementation/performancenotes.md",
            "Output quantities" => "implementation/outputquantities.md",
            "Warning: Ewald summation" => "implementation/ewaldwarning.md"
        ],
        "Usage" => [
            "Configuration files" => "usage/configfiles.md",
            "Transformation files" => "usage/transfiles.md",
            "Bloch periodicity files" => "usage/periodicfiles.md",
            "Frequency and Bloch wavevector files" => "usage/freqkfiles.md",
            "Command-line arguments" => "usage/CMDARGS.md",
            "Output file structure" => "usage/outfiles.md"
        ],
        "Examples" => [
            "Guanine + cytosine + PEC plane" => "examples/guaninecytosinePEC.md",
            "Two identical parallel graphene sheets (vacuum)" => "examples/2grapheneVAC.md"
        ],
        "API" => [
            "Command-line arguments" => "api/CMDARGS.md",
            "Exponential integral function" => "api/expintegral.md",
            "Generic matrix and vector utilities" => "api/genericmatrixvectorutils.md",
            "Maxwell Green's functions" => "api/greensfunctions.md",
            "Data structures and functions for Bloch periodicity" => "api/periodicdata.md",
            "Molecular data structures and methods" => "api/moldata.md",
            "Transformation data structures and methods" => "api/transdata.md",
            "vdW interactions and thermal radiation functions" => "api/vdWRHT.md",
        ]
    ],
    format = Documenter.HTML(prettyurls = false) # keep this line only for local builds
);
