# Frequency and Bloch wavevector files

!!! note

    Specification of Bloch wavevectors by command-line arguments or by
    input files is optional, and will be ignored when there are no
    Bloch periodic boundary conditions (discrete translational
    symmetries).

!!! note

    For vdW interaction free energies, all frequency inputs ``w`` will
    be interpreted as ``\omega = \mathrm{i}\vert w\vert``. For thermal
    radiation powers, all frequency inputs ``w`` will be interpreted
    as ``\omega = \vert w\vert``.

Frequencies and Bloch wavevectors may be supplied to the programs in
AARMBEM.jl as command-line arguments or as input files. This page goes
over the format of the latter.

## Frequency file structure

Frequencies should be specified in a file in the following structure:
```
freq1
freq2
...
freqN
```
with each value (`freq1`, `freq2`, et cetera) being replaced by a real
(currently, AARMBEM.jl does not support evaluation at arbitrary
complex frequencies, but this can be changed relatively easily in the
source code) floating-point number (possible in scientific notation),
with one value per line, and with no extraneous characters, spaces, or
blank lines. A file with a single frequency is valid. Frequencies
should be given in units of radians per second.

## Bloch wavevector file structure

Bloch wavevectors should be specified in a file in the following
structure:
```
k1x k1y k1z
k2x k2y k2z
...
kNx kNy kNz
```
with each value (`k1x`, `k1y`, et cetera) being replaced by a real
floating-point number (possible in scientific notation), with three
values per line separated by exactly a single space to represent the
Cartesian components of the Bloch wavevector, and with no extraneous
characters, spaces, or blank lines. A file with a single Bloch
wavevector is valid, but it must contain all three components
explicitly. Wavevectors should be given in units of radians per meter.