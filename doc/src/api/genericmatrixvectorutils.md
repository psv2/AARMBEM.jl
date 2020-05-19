# Generic matrix and vector utilities

These functions can all be used independently of AARMBEM.jl.

## Involving a full matrix and a diagonal matrix as a vector

```@docs
addtodiag!
addtodiag
subtractfromdiag!
subtractfromdiag
mulMatDiagleft!
mulMatDiagleft
mulMatDiagright!
mulMatDiagright
```

## Involving only a vector

```@docs
vecNto3N
```

## Involving only a square matrix

```@docs
sympart!
sympart
hermpart!
hermpart
ahermpart!
ahermpart
```