# Molecular configuration files

A major ingredient of running the command-line programs in AARMBEM.jl
is the molecular configuration file. The structure of such a file is
described in the following sections. As discussed in the page about
command-line arguments, these calculations can be done including the
effects of phonons (arising from displacements of the screened nuclear
oscillators), but they can also be done in the absence of phonons,
effectively taking the valence electronic oscillators to be coupled to
fixed, rigid nuclei.

## Configuration file structure

The configuration file may have any file type, but it must be readable
as plain text in a typical text editor. Given this, it is typical for
the configuration file to have the `.txt` file type.

For ``N_{\mathrm{mol}}`` molecular bodies, the configuration file is
made of ``N_{\mathrm{mol}}`` blocks defining the properties of each
molecular body. A block is structured as in the following example.

```
MOL mol1
Qe Qefile.txt
Me Mefile.txt
Ke Kefile.txt
Be Befile.txt
MI MIfile.txt
BI BIfile.txt
KI KIfile.txt
xyz xyzfile.txt
Bescale 1.0
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL
```

In terms of formatting of this text block, leading and trailing spaces
as well as blank lines are ignored. However, there are certain
features that cannot be ignored or changed.

The first line defining a block must be the word `MOL` (case
sensitive), followed by a space, followed by a string defining the
molecular label (in this example `mol1`); this label string can be
anything, but should, for ease of use, be somehow descriptive of the
molecular body. The last line, before any new block starting with
`MOL` is defined, must be `ENDMOL`.

The lines in between can be put in any order. For this molecule to be
labeled ``n``, the file used for constructing the matrices
``Q_{\mathrm{e}n}`` is specified as `Qefile.txt`, coming in the same
line as `Qe` followed by exactly one space; no nonwhitespace
characters are present after the specification of the file name. The
corresponding files for ``M_{\mathrm{e}n}``, ``K_{\mathrm{e}n}``,
``B_{\mathrm{e}n}``, ``M_{\mathrm{I}n}``, ``B_{\mathrm{I}n}``,
``K_{\mathrm{I}n}``, and the atomic coordinates come respectively as
`Mefile.txt`, `Kefile.txt`, `Befile.txt`, `MIfile.txt`, `BIfile.txt`,
`KIfile.txt`, and `xyzfile.txt`, each respectively following the
keywords `Me`, `Ke`, `Be`, `MI`, `BI`, `KI`, or `xyz` and (in each
case) following exactly one space.

Of the above keywords, with or without phonons, the lines with `Qe`,
`Me`, `Ke`, and `xyz` are required. With phonons, the keywords
corresponding to the lines with `MI` and `KI` are required. With or
without phonons, the keywords corresponding to the lines with `Be` and
`BI` are always optional. If the line with `Be` is omitted,
``B_{\mathrm{e}n}`` is set to vanish, and likewise in the case with
phonons, if the line with `BI` is omitted, ``B_{\mathrm{I}n}`` is set
to vanish; that said, if both of these vanish, thermal radiation
powers will vanish, as those require nontrivial material dissipation.

The line `Bescale bescaleval` allows for easily multiplying the matrix
``B_{\mathrm{e}n}`` by a uniform numerical factor `bescaleval`, which
is a floating-point number that must correspond to a real nonnegative
value (and may be written in scientific notation): that value
`bescaleval` must follow the keyword `Bescale` and then exactly one
space. In the above example, `bescaleval` is `1.0`. This is an
optional keyword, and in its absence, ``B_{\mathrm{e}n}`` is unchanged
from what is specified through the file associated with the keyword
`Be`.

In the case with phonons, the line beginning with the keyword
`KIconst` specifies an artificial spring constant to constrain one or
two screened nuclear oscillators. This is useful for the following
reason, as detailed later in this page. A physical spring constant
matrix ``K_{\mathrm{I}n}`` for a collection of ``N_{n}`` atoms in 3
spatial dimensions will be real-symmetric positive-semidefinite, with
exactly 6 eigenvalues that are zero, with 3 modes corresponding to
inertial motion of the center of mass in each of the 3 Cartesian
directions, and 3 modes corresponding to rigid rotation via Euler
angles about the center of mass (though for an atomically-thin linear
molecular body, like a carbyne wire, one of the rotational modes
disappears, so there would only be 5 eigenvalues that vanish). These
zero-frequency eigenvalues lead to divergence of the polarizability
matrix ``\alpha_{n}`` as ``\omega \to 0`` (along the real or imaginary
axes), which could potentially lead to problems particularly for vdW
interaction free energies where the integrand at ``\omega = 0`` does
not vanish. A useful way of squeezing out these modes is to
artificially constrain two (for compact molecular bodies, or only one
for molecular bodies with commensurate discrete translational
symmetry) screened nuclear oscillators with extra spring constants;
for compact molecular bodies, constraining the screened nuclear
oscillators that are farthest apart from each other can effectively
preserve the phonon dispersion properties of most of the molecule, but
for periodically extended molecular bodies, constraints are generally
not recommended as their periodic repetition for all unit cells will
lead to spurious modifications to the material properties.

The artificial spring constant constraints may be specified in a few
different ways; in each case, the value is given in units of newtons
per meter. Also, this is an optional keyword: if that line is absent,
the artificial spring constant is set by default to vanish. Following
the keyword `KIconst` and then exactly one space, exactly 1, 3, 6, or
9 values may be specified. If the line is specified as `KIconst
KIconst1`, where `KIconst1` is a real nonnegative number specified as
a floating-point value (possibly with scientific notation), then the
spring constant tensor is simply `KIconst1` multiplied by the 3-by-3
identity tensor, implying an isotropic coupling; in the example above,
`KIconst1` is `1.0e3` (i.e. ``10^{3}~\mathrm{N/m}``). If that value
happens to be negative, its absolute value is taken. If the line is
specified as `KIconst KIconst1 KIconst2 KIconst3`, then the values
`KIconst1`, `KIconst2`, and `KIconst3` are respectively assigned as
the ``xx``-, ``yy``-, and ``zz``-components of a tensor that is
diagonal in the basis of the Cartesian axis unit vectors; if any of
these values are negative, their absolute values are used. If the line
is specified as `KIconst KIconst1 KIconst2 KIconst3 KIconst4 KIconst5
KIconst6`, then the values `KIconst1`, `KIconst2`, `KIconst3`,
`KIconst4`, `KIconst5`, and `KIconst6` are respectively assigned as
the ``xx``-, ``xy``-, ``xz``-, ``yy``-, ``yz``-, and ``zz``-components
of a general tensor, with the ``yx``-, ``zx``-, and ``zy``-components
respectively assigned to be the same as the ``xy``-, ``xz``-, and
``yz``-components to ensure that the tensor is real-symmetric; once
this is done, the eigenvalues and eigenvectors are found, any
eigenvalues that are negative are replaced by their absolute values,
and these are then multiplied by the corresponding eigenvector
projection operators and added again to yield a real-symmetric
positive-definite tensor spring constant constraint along general
orthonormal axes in 3 dimensions. If the line is specified as `KIconst
KIconst1 KIconst2 KIconst3 KIconst4 KIconst5 KIconst6 KIconst7
KIconst8 KIconst9`, then the values `KIconst1`, `KIconst2`,
`KIconst3`, `KIconst4`, `KIconst5`, `KIconst6`, `KIconst7`,
`KIconst8`, and `KIconst9` are respectively assigned as the ``xx``-,
``xy``-, ``xz``-, ``yx``- ``yy``-, ``yz``-, ``zx``-, ``zy``-, and
``zz``-components of a general tensor, which is then symmetrized,
following which (as in the case with 6 values specified) its
eigenvectors are preserved but its eigenvalues are replaced by their
absolute values (if needed) to yield a real-symmetric
positive-definite tensor spring constant constraint along general
orthonormal axes in 3 dimensions.

The lines `KIconstatom1 KIconstatom1val` and `KIconstatom2
KIconstatom2val` specify which atoms are to be constrained. Both
`KIconstatom1val` and `KIconstatom2val` are to be specified as
integers between 1 and ``N_{n}`` if the molecular body labeled ``n``
has ``N_{n}`` atoms (in the example above, `KIconstatom1val` is `1`
and `KIconstatom2val` is `2`): any specification outside of this
range, or any issues like `KIconstatom1val` being larger than or equal
to `KIconstatom2val`, are dealt with by appropriately setting
`KIconstatom1val` to 1 or `KIconstatom2val` to ``N_{n}``. These are
optional keywords: in their absence, the aforementioned default values
of 1 and ``N_{n}`` are used. Furthermore, for Bloch periodic boundary
conditions, only `KIconstatom1val` is used. Currently, AARMBEM.jl uses
the same artificial spring constant constraint for both atoms in a
molecular body, and does not support setting separate constraints for
different screened nuclear oscillators in a given molecular body.

The only time when a molecular body may be specified without all of
this data is if it is a duplicate of another molecular body. Such a
specification may be written as
```
MOL mol2
DUPLICATE mol1
ENDMOL
```
where the keyword ``DUPLICATE`` (case sensitive) is followed by a
space and then the exact string label of another molecular body, to
say that these share all of their material properties and only differ
in position or orientation. Any number of molecular bodies may be
listed as duplicates of a given molecule body. However, the molecular
label string must not itself be a duplicate (i.e. it must correspond
to a molecular body whose properties are explicitly specified), and
molecular bodies with explicitly specified properties must be listed
before molecular bodies that are duplicates of those. For instance,
the following is allowed:
```
MOL mol1
Qe Qefile1.txt
Me Mefile1.txt
Ke Kefile1.txt
Be Befile1.txt
MI MIfile1.txt
BI BIfile1.txt
KI KIfile1.txt
xyz xyzfile1.txt
Bescale 1.0
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL

MOL mol2
DUPLICATE mol1
ENDMOL

MOL mol3
DUPLICATE mol1
ENDMOL

MOL mol4
Qe Qefile4.txt
Me Mefile4.txt
Ke Kefile4.txt
Be Befile4.txt
MI MIfile4.txt
BI BIfile4.txt
KI KIfile4.txt
xyz xyzfile4.txt
Bescale 1.0e4
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL

MOL mol5
DUPLICATE mol1
ENDMOL

MOL mol6
DUPLICATE mol4
ENDMOL

MOL mol7
DUPLICATE mol1
ENDMOL
```
but neither the following:
```
MOL mol1
Qe Qefile1.txt
Me Mefile1.txt
Ke Kefile1.txt
Be Befile1.txt
MI MIfile1.txt
BI BIfile1.txt
KI KIfile1.txt
xyz xyzfile1.txt
Bescale 1.0
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL

MOL mol2
DUPLICATE mol1
ENDMOL

MOL mol3
DUPLICATE mol4
ENDMOL

MOL mol4
Qe Qefile4.txt
Me Mefile4.txt
Ke Kefile4.txt
Be Befile4.txt
MI MIfile4.txt
BI BIfile4.txt
KI KIfile4.txt
xyz xyzfile4.txt
Bescale 1.0e4
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL

MOL mol5
DUPLICATE mol1
ENDMOL

MOL mol6
DUPLICATE mol4
ENDMOL

MOL mol7
DUPLICATE mol1
ENDMOL
```
nor the following:
```
MOL mol1
Qe Qefile1.txt
Me Mefile1.txt
Ke Kefile1.txt
Be Befile1.txt
MI MIfile1.txt
BI BIfile1.txt
KI KIfile1.txt
xyz xyzfile1.txt
Bescale 1.0
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL

MOL mol2
DUPLICATE mol1
ENDMOL

MOL mol3
DUPLICATE mol2
ENDMOL

MOL mol4
Qe Qefile4.txt
Me Mefile4.txt
Ke Kefile4.txt
Be Befile4.txt
MI MIfile4.txt
BI BIfile4.txt
KI KIfile4.txt
xyz xyzfile4.txt
Bescale 1.0e4
KIconst 1.0e3
KIconstatom1 1
KIconstatom2 2
ENDMOL

MOL mol5
DUPLICATE mol1
ENDMOL

MOL mol6
DUPLICATE mol4
ENDMOL

MOL mol7
DUPLICATE mol1
ENDMOL
```
are allowed.

## Input file structure

Here, the structures of the input files following the keywords `Qe`,
`Me`, `Ke`, and so on are described.

The files following the keywords `Qe`, `Me`, and so on can have any
file type, but should be readable as plain text. In particular, the
files following the keywords `Qe`, `Me`, `Ke`, and if present, `Be`,
`MI`, and `BI` must a list of ``N_{n}`` values if molecular body ``n``
has ``N_{n}`` atoms, formatted as follows:
```
val1
val2
val3
...
valNn
```
so each value should be specified as a real nonnegative floating-point
quantity (possibly in scientific notation), and there should be no
spaces, blank lines, or extraneous characters. The values in the files
following the keywords `Me` and `MI` should be in units of kilograms,
those in the files following the keywords `Be` and `BI` should be in
units of kilograms per second, those in the file following the keyword
`Ke` should be in units of newtons per meter, and those in the file
following the keyword `Qe` should be in
``\mathrm{kg}^{1/2} \cdot \mathrm{m}^{3/2} \cdot \mathrm{s}^{-1}``
obtained by dividing the charge in coulombs by the square root of the
vacuum permittivity in SI.

The file following the keyword `xyz` can have any file type, but
should be readable as plain text. In particular, the file following
the keyword `xyz` must be a list of 3 lines, each of which has
``N_{n}`` values separated only by commas, if molecular body ``n`` has
``N_{n}`` atoms, formatted as follows:
```
r_{1x},r_{2x},...,r_{Nnx}
r_{1y},r_{2y},...,r_{Nny}
r_{1z},r_{2z},...,r_{Nnz}
```
so each value ``r_{pi}`` should be specified as a real floating-point
quantity (possibly in scientific notation), and there should be no
spaces, blank lines, or extraneous characters. The values in the file
following the keyword `xyz` should be in units of meters.

The file following the keyword `KI` can have any file type, but should
be readable as plain text. In particular, the file following the
keyword `xyz` for a compact molecular body labeled ``n`` must be a
list of ``3N_{n}`` lines, each of which has ``3N_{n}`` values
separated only by commas, if molecular body ``n`` has ``N_{n}`` atoms,
formatted as follows:
```
K_{I1x,1x},K_{I1x,1y},K_{I1x,1z},K_{I1x,2x},K_{I1x,2y},K_{I1x,2z},...,K_{I1x,Nnx},K_{I1x,Nny},K_{I1x,Nnz}
K_{I1y,1x},K_{I1y,1y},K_{I1y,1z},K_{I1y,2x},K_{I1y,2y},K_{I1y,2z},...,K_{I1y,Nnx},K_{I1y,Nny},K_{I1y,Nnz}
K_{I1z,1x},K_{I1z,1y},K_{I1z,1z},K_{I1z,2x},K_{I1z,2y},K_{I1z,2z},...,K_{I1z,Nnx},K_{I1z,Nny},K_{I1z,Nnz}
K_{I2x,1x},K_{I2x,1y},K_{I2x,1z},K_{I2x,2x},K_{I2x,2y},K_{I2x,2z},...,K_{I2x,Nnx},K_{I2x,Nny},K_{I2x,Nnz}
K_{I2y,1x},K_{I2y,1y},K_{I2y,1z},K_{I2y,2x},K_{I2y,2y},K_{I2y,2z},...,K_{I2y,Nnx},K_{I2y,Nny},K_{I2y,Nnz}
K_{I2z,1x},K_{I2z,1y},K_{I2z,1z},K_{I2z,2x},K_{I2z,2y},K_{I2z,2z},...,K_{I2z,Nnx},K_{I2z,Nny},K_{I2z,Nnz}
...,...,...,...,...,...,...,...,...,...
K_{INnx,1x},K_{INnx,1y},K_{INnx,1z},K_{INnx,2x},K_{INnx,2y},K_{INnx,2z},...,K_{INnx,Nnx},K_{INnx,Nny},K_{INnx,Nnz}
K_{INny,1x},K_{INny,1y},K_{INny,1z},K_{INny,2x},K_{INny,2y},K_{INny,2z},...,K_{INny,Nnx},K_{INny,Nny},K_{INny,Nnz}
K_{INnz,1x},K_{INnz,1y},K_{INnz,1z},K_{INnz,2x},K_{INnz,2y},K_{INnz,2z},...,K_{INnz,Nnx},K_{INnz,Nny},K_{INnz,Nnz}
```
so each value ``K_{Ipi,qj}`` should be specified as a real
floating-point quantity (possibly in scientific notation), and there
should be no spaces, blank lines, or extraneous characters. The values
in the file following the keyword `KI` should be in units of newtons
per meter.

There are a few points of note for the file following the keyword
`KI`, if present. For a compact molecular body, it must correspond to
a real-symmetric positive-semidefinite matrix, with exactly 6
eigenvalues (or 5 for a purely linear 1-dimensional atomically thin
wire-like molecular body) that are effectively zero (i.e. are
significantly smaller in magnitude than the remaining eigenvalues,
though this is left to the judgment of the user). For a molecular body
obeying periodic boundary conditions, there may or may not be nonzero
harmonic couplings between screened nuclear oscillators across
different unit cells. If there are not, then the above format for the
file following the keyword `KI` should be used. If there are, then a
different format should be used. In particular, for 1 periodic
dimension, given a unit cell (which may not be the primitive unit
cell) with primitive lattice vector ``\vec{a}``, any lattice vector
may be represented as ``\vec{R} = n\vec{a}`` for integer ``n``
(separate from the molecular labels). Likewise, for 2 periodic
dimensions, given a unit cell (which may not be the primitive unit
cell) with primitive lattice vectors ``\vec{a}_{1}`` and
``\vec{a}_{2}`` (which might not be orthogonal), any lattice vector
may be represented as ``\vec{R} = n_{1} \vec{a}_{1} + n_{2}
\vec{a}_{2}`` for integers ``n_{1}, n_{2}``. In either case, given a
base unit cell labeled 0 corresponding to the lattice vector ``\vec{R}
= 0``, the matrix block ``K_{\mathrm{I}(\vec{R}, 0)}`` describes
couplings from screened nuclear oscillators in unit cell 0 to those in
unit cell ``\vec{R}`` (i.e. ``(K_{\mathrm{I}(\vec{R}, 0)})_{pi,qj}``
describes the coupling of the screened nuclear oscillator in atom
``q`` in unit cell 0 along Cartesian direction ``j`` to the periodic
image of the screened nuclear oscillator in atom ``p`` in unit cell
``\vec{R}`` along Cartesian direction ``i``). Because the unit cell
indexing is always of the form ``(\vec{R}, 0)`` (i.e. atom ``q``
always lies in unit cell 0), notation may be slightly abused such that
in 1 periodic dimension, ``K_{\mathrm{I}(n)}`` is the matrix block
``K_{\mathrm{I}(n\vec{a}, 0)}``, and in 2 periodic dimensions,
``K_{\mathrm{I}(n_{1}, n_{2})}`` is the matrix block
``K_{\mathrm{I}(n_{1} \vec{a}_{1} + n_{2} \vec{a}_{2}, 0)}``. Each of
the blocks ``K_{\mathrm{I}(n)}`` or ``K_{\mathrm{I}(n_{1}, n_{2})}``
should be ``3N_{n} \times 3N_{n}`` (where ``n`` in the context of
``N_{n}`` is again the molecular label) matrices. Thus, for
internuclear couplings in either 1 periodic dimension or 2 periodic
dimensions where couplings are only nontrivial along ``\vec{a}_{1}``,
the matrix blocks should be arranged (with no blank lines in the file)
as
```math
\begin{bmatrix}
K_{\mathrm{I}(-n_{1\max}, 0)} \\
\vdots \\
K_{\mathrm{I}(-1, 0)} \\
K_{\mathrm{I}(0, 0)} \\
K_{\mathrm{I}(1, 0)} \\
\vdots \\
K_{\mathrm{I}(n_{1\max}, 0)} 
\end{bmatrix}
```
(replacing ``K_{\mathrm{I}(n_{1}, 0)}`` in 2 periodic dimensions with
``K_{\mathrm{I}(n)}`` in 1 periodic dimension) where ``n_{1\max}`` (or
``n_{\max}``) is the largest unit cell index in either direction along
that periodic dimension where internuclear couplings between any atom
``p`` in that unit cell and any atom ``q`` in unit cell 0 are
nonzero. In 2 periodic dimensions where couplings are only nontrivial
along ``\vec{a}_{2}``, the matrix blocks should be arranged (through
horizontal concatenation with commas as usual, with no other
extraneous characters, spaces, or blank lines) as
```math
\begin{bmatrix}
K_{\mathrm{I}(0, -n_{2\max})} & \ldots & K_{\mathrm{I}(0, -1)} & K_{\mathrm{I}(0, 0)} & K_{\mathrm{I}(0, 1)} & \ldots & K_{\mathrm{I}(0, n_{2\max})} 
\end{bmatrix}
```
where ``n_{2\max}`` is the largest unit cell index in either direction
along that periodic dimension where internuclear couplings between any
atom ``p`` in that unit cell and any atom ``q`` in unit cell 0 are
nonzero. In 2 periodic dimensions where couplings are nontrivial along
both primitive lattice vector directions ``\vec{a}_{1}`` and
``\vec{a}_{2}``, the matrix blocks should be arranged (through
horizontal concatenation with commas as usual and direct vertical
concatenation, with no other extraneous characters, spaces, or blank
lines) as
```math
\begin{bmatrix}
K_{\mathrm{I}(-n_{1\max}, -n_{2\max})} & \ldots & K_{\mathrm{I}(-n_{1\max}, -1)} & K_{\mathrm{I}(-n_{1\max}, 0)} & K_{\mathrm{I}(-n_{1\max}, 1)} & \ldots & K_{\mathrm{I}(-n_{1\max}, n_{2\max})} \\
\vdots & \ddots & \vdots & \vdots & \vdots & \ddots & \vdots \\
K_{\mathrm{I}(-1, -n_{2\max})} & \ldots & K_{\mathrm{I}(-1, -1)} & K_{\mathrm{I}(-1, 0)} & K_{\mathrm{I}(-1, 1)} & \ldots & K_{\mathrm{I}(-1, n_{2\max})} \\
K_{\mathrm{I}(0, -n_{2\max})} & \ldots & K_{\mathrm{I}(0, -1)} & K_{\mathrm{I}(0, 0)} & K_{\mathrm{I}(0, 1)} & \ldots & K_{\mathrm{I}(0, n_{2\max})} \\
K_{\mathrm{I}(1, -n_{2\max})} & \ldots & K_{\mathrm{I}(1, -1)} & K_{\mathrm{I}(1, 0)} & K_{\mathrm{I}(1, 1)} & \ldots & K_{\mathrm{I}(1, n_{2\max})} \\
\vdots & \ddots & \vdots & \vdots & \vdots & \ddots & \vdots \\
K_{\mathrm{I}(n_{1\max}, -n_{2\max})} & \ldots & K_{\mathrm{I}(n_{1\max}, -1)} & K_{\mathrm{I}(n_{1\max}, 0)} & K_{\mathrm{I}(n_{1\max}, 1)} & \ldots & K_{\mathrm{I}(n_{1\max}, n_{2\max})} 
\end{bmatrix}
```
where ``n_{1\max}`` and ``n_{2\max}`` are the largest unit cell
indices in either direction along that corresponding periodic
dimension where internuclear couplings between any atom ``p`` in that
unit cell and any atom ``q`` in unit cell 0 are nonzero.

Whether for compact or Bloch periodic molecular bodies, the matrix
``K_{\mathrm{I}}`` in real space must also satisfy two other
properties. Every off-diagonal ``3\times 3`` tensor
``\mathbb{K}_{\mathrm{I}pq}`` for atoms ``p \neq q`` (and for Bloch
periodic molecular bodies, this includes tensor blocks in the block
matrices ``K_{\mathrm{I}(\vec{R}, 0)}`` for ``\vec{R} \neq 0``, so
atom ``q`` would be in unit cell 0 while atom ``p`` may be in another
unit cell ``\vec{R}``) must be real-symmetric
negative-definite. Furthermore, the ``3\times 3`` tensors in the
diagonal blocks must satisfy ``K_{\mathrm{I}qi,qj} = -\sum_{p}
K_{\mathrm{I}pi,qj}``, (for Cartesian tensor indices ``i, j``) where
for a compact molecular body atoms ``p, q`` are within the body itself
while for a periodic molecular body with nontrivial couplings to other
unit cells ``q`` lies within unit cell 0 while ``p`` runs over all
other atoms in all unit cells (including the periodic image of ``q``
in other unit cells).

!!! note
    
    For all of these data files, AARMBEM.jl does not check that the
    corresponding lists have the correct number of elements, and will
    simply throw an error if the formatting or list lengths are
    incorrect (i.e. it does not try to correct the error).

!!! warning
    
    As a reminder, AARMBEM.jl does not perform the DFT calculations
    needed to generate the relevant molecular input data files. This
    must be done separately.

!!! warning
    
    AARMBEM.jl does not check whether the matrix ``K_{\mathrm{I}}`` in
    real space obeys the rules of positive-definiteness overall,
    negative-definiteness for the off-diagonal ``3\times 3`` tensor
    blocks, or the summation constraint for the diagonal ``3\times 3``
    tensor blocks. These need to be verified, and enforced if need be,
    beforehand. Furthermore, if ``K_{\mathrm{I}}`` is computed via
    DFT, such calculations will typically yield nontrivial couplings
    between all pairs of atoms at arbitrary distances (though the
    magnitudes will decrease significantly with increasing distance),
    and the results may lead to numerical error depending on the
    extent to which any of these conditions are violated; as a result,
    it may be necessary to symmetrize off-diagonal ``3\times 3``
    tensor blocks manually, set couplings beyond a certain bond
    distance to zero, and recompute the ``3\times 3`` diagonal tensor
    blocks manually.