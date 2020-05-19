# Molecular Bloch periodicity files

An ingredient of running the command-line programs in AARMBEM.jl is
the file specifying Bloch periodicity. This is optional, and should
only be used if Bloch periodicity is present. If this file is not
supplied, the default behavior is to assume no Bloch periodicity
(i.e. compact molecular bodies).

The Bloch periodicity file may have any file type, but it must be
readable as plain text in a typical text editor. Given this, it is
typical for the Bloch periodicity file to have the `.txt` file type.

If there is discrete translational symmetry in only 1 dimension, the
file should simply have 3 values
```
ax
ay
az
```
where the values `ax`, `ay`, and `az` are replaced by real
floating-point numbers (possibly in scientific notation) representing
the Cartesian components of the lattice vector ``\vec{a}``.

If there is discrete translational symmetry in 1 dimensions, the file
should have 6 values
```
a1x
a1y
a1z
a2x
a2y
a2z
```
where the values `a1x`, `a1y`, `a1z`, `a2x`, `a2y`, and `a2z` are
replaced by real floating-point numbers (possibly in scientific
notation) representing the Cartesian components of the respective
lattice vector ``\vec{a}_{1}`` and ``\vec{a}_{2}``.

In both cases, there should be no extraneous characters, spaces, or
blank lines.