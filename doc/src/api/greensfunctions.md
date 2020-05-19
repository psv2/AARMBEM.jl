# Maxwell Green's functions

```@docs
c
```

## Vacuum Green's function

This section of functions can be used outside of AARMBEM.jl.

```@docs
GFVACGG!
GFVACGG
GFVACGG_sca!
GFVACGG_sca
GFVACGG0!
GFVACGG0
GFVACGG0_sca!
GFVACGG0_sca
GFVACGGcoincident!
GFVACGGcoincident
GFVACGG1Ewald!
GFVACGG1Ewald
GFVACGG2Ewald!
GFVACGG2Ewald
GFVACGGEwaldSR
GFVACGGEwaldLR1D
GFVACGGEwaldLR2D
```

## Scattering Green's function above a PEC plane at ``z = 0``

This section of functions can be used outside of AARMBEM.jl.

```@docs
GFPECGG!
GFPECGG
GFPECGG0!
GFPECGG0
GFPECGG1Ewald!
GFPECGG1Ewald
GFPECGG2Ewald!
GFPECGG2Ewald
```

## Green's function matrix assembly routines

```@docs
assemble1MolGvacinf!
assemble1MolGvacinfdiag!
subtract1MolalphainvGvac!
assembleGscaMolBlock!
assembleGvacMolBlock!
```