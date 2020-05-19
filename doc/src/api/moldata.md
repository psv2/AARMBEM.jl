# Molecular data structures and methods

## Data structures

```@docs
OneMol
OneMol_WithPhonons
OneMol_NoPhonons
MolSystem
```

## Frequency-independent data structure initialization functions

```@docs
read1MolDataWithPhonons
read1MolDataNoPhonons
readMolSystem
```

## Frequency-dependent polarizability calculation functions

```@docs
constructMolAlpha!(myMolData::OneMol_WithPhonons{FT}, myPeriodicData::PeriodicData{FT}, freq::Union{FT, Complex{FT}}, k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat
constructMolAlpha!(myMolData::OneMol_NoPhonons{FT}, myPeriodicData::PeriodicData{FT}, freq::Union{FT, Complex{FT}}, k::Array{FT, 1}=zeros(FT, 3)) where FT<:AbstractFloat
KIk
```