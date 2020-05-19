# Molecular transformation files

A major ingredient of running the command-line programs in AARMBEM.jl
is the molecular transformation file. For computation of vdW
interaction free energies and thermal radiation powers, it is
frequently useful to translate and rotate different molecular bodies
in different orientations, whether in vacuum or above a PEC plane at
``z = 0``. The following sections explain the different parts of the
transformation file and their structure. The use of this file is
technically optional, in that if this file is not present, then all
molecular bodies are set to have their unweighted centers of mass at
the origin and orientations unchanged from those respectively
specified in the files associated with the keyword `xyz` in the
molecular configuration file. However, this would lead to multiple
molecular bodies overlapping with each other or with a PEC plane (if
present). Therefore, this file should be used in almost all cases.

The transformation file may have any file type, but it must be
readable as plain text in a typical text editor. Given this, it is
typical for the transformation file to have the `.txt` file type.

## Baseline positions and orientations

The first section of the file specifies the baseline unweighted center
of mass positions and orientations of all molecular bodies. As a
generic example, it is structured as
```
BASELINE
mol1 x1 y1 z1
mol2 x2 y2 z2 phi2 theta2 psi2
mol4 x4 y4 z4 phi4 theta4 psi4
mol5 x5 y5 z5
ENDBASELINE
```
and must have the first line be `BASELINE` and the last line be
`ENDBASELINE` (both case sensitive with no other spaces). In each
line, the first string is the molecular label string for the
corresponding molecular body (as specified in the molecular
configuration file: the label string must exactly match). Following
that and exactly one space can be 3 or 6 values (each separated by
exactly one space, with no extraneous characters). The first, second,
and third of those values respectively correspond to the ``x``-,
``y``-, and ``z``- coordinates of the unweighted center of mass. The
fourth, fifth, and sixth of those values respectively correspond to
the Euler angles ``\varphi``, ``\theta``, and ``\psi`` in *degrees*
relative to the orientation specified in the molecular input data file
for atomic coordinates. All of these are real floating point values
(possibly in scientific notation). If only three values are present,
then the orientation is unchanged from that specified in the molecular
input data file for atomic coordinates, and the three values present
correspond to the unweighted center of mass. If a molecular body is
not specified in a line in this section, it is assumed to have a
default unweighted center of mass at the orientation and a default
orientation unchanged from that specified in the molecular input data
file for atomic coordinates.

## Transformations

The second section of the file specifies the baseline unweighted
center of mass positions and orientations of all molecular bodies
*relative to their baseline centers of mass and orientations*. As a
generic example, a transformation is structured as
```
TRANS mytrans
mol1 x1 y1 z1
mol4 x4 y4 z4 phi4 theta4 psi4
ENDTRANS
```
with the first line being the keyword `TRANS` (case sensitive)
followed by exactly one space and then a string labeling the
transformation (in this case `mytrans`: it should usually be
descriptive or enumerated in some way for ease of extracting data
later) and the last line being the keyword `ENDTRANS` (case
sensitive). The lines in between are similar to the specification of
the baseline centers of mass and orientations of each molecular body,
again emphasizing the caveat that the positions and angles are
actually Cartesian and angular *displacements* from the baseline
values.