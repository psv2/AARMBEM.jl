# Data structures and rigid transformations

This page describes the data structures used in AARMBEM.jl, as well as
conventions for rigid transformations of molecular bodies. The basic
idea is that for many EM scattering calculations, including vdW
interaction energies and thermal radiation powers, data is needed at
each frequency (and Bloch wavevector, if appropriate) and for many
different geometric configurations, in which each geometric
configuration consists of a set of rigid transformations applied to
molecular and macroscopic bodies; as AARMBEM.jl currently only
implements vacuum or the presence of a PEC plane at ``z = 0`` as the
macroscopic body, all rigid transformations refer to those applied to
molecular bodies. The data for each molecular body in isolation is
independent of the rigid transformation, and therefore may be computed
once per frequency (and Bloch wavevector, if appropriate); this has
the additional benefit of being able to reuse data if multiple
molecular bodies are used which are copies of each other in
isolation. After this, data is computed for the Green's function and
T-operator matrices ``G^{\mathrm{mac}}`` and ``T`` for each rigid
transformation (which may include translations as well as proper
rotations), efficiently reusing contributions from molecular data
computed in isolation when possible. Further details are given in the
remainder of this page.

## Molecular data in isolation

In the RMB framework, every molecular body ``n`` is characterized by
properties in isolation, including the number of atoms ``N_{n}``,
atomic coordinates ``\{\vec{r}_{p}\}``, oscillator parameters, and
atomic polarizabilities. In particular, as the charges, masses,
damping coefficients, and valence electronic oscillator spring
constants are all isotropic, they can efficiently be stored as
``N_{n}``-element vectors; only the screened nuclear oscillator spring
constants must be stored as a full ``3N_{n} \times 3N_{n}``
matrix. From these quantities, the polarizability matrix
``\alpha_{n}`` may be computed once per frequency (and Bloch
wavevector, if appropriate). Furthermore, the diagonal molecular block
of the vacuum Green's function matrix ``G^{(0)}_{nn}`` for each
molecular body ``n`` likewise only depends on the relative positions
of atoms in that body, not on absolute positions in space, and is
therefore independent of rigid transformations, requiring computation
only once per frequency (and Bloch wavevector, if appropriate).

Each molecular body can also be specified to have a certain
*unweighted* center of mass position and a certain spatial orientation
relative to the input data. Frequently, it may be desired to perform
computations involving multiple molecular bodies that have entirely
identical properties. It would be wasteful to duplicate molecular data
in such a case, so the AARMBEM.jl code stores only data for unique
molecular bodies, while separately storing only the centers of mass
and orientations for each molecular body (distinguishing between
otherwise identical molecular bodies). The polarizability and vacuum
Green's function matrices ``\alpha_{n}`` and ``G^{(0)}_{nn}`` for each
molecular body ``n`` in isolation are reused to construct the Green's
function and T-operator matrices ``G^{\mathrm{mac}}`` and ``T``, using
appropriate rotation matrices to reproduce the desired orientations
while translating the atomic coordinates back and forth accordingly.

## Rigid transformations

Each molecular body has a baseline center of mass and a baseline
orientation; if none are specified, the baseline center of mass is
chosen to be the origin, and the baseline orientation is chosen to be
unchanged from that corresponding to the input data for the atomic
coordinates for that body. Each molecular body may also undergo
translations of the center of mass and rigid rotations about the
center of mass; the translation vectors are relative to the baseline
center of mass, and the rotations are relative to the baseline
orientation.

Rotations are specified through Euler angles ``(\varphi, \theta,
\psi)``. The particular convention used by AARMBEM.jl is that an
atomic coordinate ``\vec{r}``, represented by the column vector
``\begin{bmatrix} x, y, z \end{bmatrix}^{\top}``, is transformed by a
rotation into ``R_{z} (-\psi) R_{y} (-\theta) R_{x} (-\varphi)
\begin{bmatrix} x, y, z \end{bmatrix}^{\top}``. The rotation matrices
are defined as
```math
R_{x} (-\varphi) = \begin{bmatrix}
1 & 0 & 0 \\
0 & \cos(\varphi) & \sin(\varphi) \\
0 & -\sin(\varphi) & \cos(\varphi)
\end{bmatrix} \\
R_{y} (-\theta) = \begin{bmatrix}
\cos(\theta) & 0 & -\sin(\theta) \\
0 & 1 & 0 \\
\sin(\theta) & 0 & \cos(\theta)
\end{bmatrix} \\
R_{z} (-\psi) = \begin{bmatrix}
\cos(\psi) & \sin(\psi) & 0 \\
-\sin(\psi) & \cos(\psi) & 0 \\
0 & 0 & 1
\end{bmatrix}
```
in terms of these Euler angles.

The individual molecular matrices ``\alpha_{n}`` and ``G^{(0)}_{nn}``
are appropriately rotated before stamping into the diagonal blocks of
``G^{\mathrm{mac}}`` or ``T``. This also accounts for the possibility
that there may be molecular bodies with identical material properties
but with different centers of mass or orientations (with only the
latter having any effect on these matrices, with the correction being
given through application of these rotation matrices as
appropriate). This allows for efficient reuse of data, as computing
the Green's function matrix elements is typically the computational
bottleneck. Additionally, if a set of molecular bodies is not affected
from one rigid transformation to the next, their diagonal and
off-diagonal matrix blocks are not changed either. Rigid
transformations are performed in a sorted order (regardless of the
order specified), such that the displacement of the last molecular
body changes first, then the displacement of the next-to-last
molecular body, and so on, and then the orientation of the last
molecular body, then the orientation of the next-to-last molecular
body, and so on; this means that if there are many changes in position
but few changes in orientation, as is often the case in many of these
calculations, fewer transformations need to be applied and
consequently fewer matrix blocks need to be recomputed.