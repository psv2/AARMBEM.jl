# Theoretical background

This page explains the theoretical background for the RMB framework,
and for its implementation in AARMBEM.jl. After explaining notation,
the main formulas and assumptions are summarized, and then the
implementation in AARMBEM.jl is explained; more details may be found
in work by Venkataram et al [^VenkataramARXIV2020].

## Notation

The following is a table of common notation used throughout this site,
and particularly this page.

| Notation | Meaning | Comments |
|:---------|:--------|:---------|
| ``\vert\vec{v}\rangle`` | ``v_{i} (\vec{x})`` | semiclassical field, *not* second-quantized state |
| ``\langle \vec{u}, \vec{v} \rangle`` | ``\sum_{i} \int u_{i}^{\star} (\vec{x}) v_{i} (\vec{x})~\mathrm{d}^{3} x`` | basis-independent conjugated inner product |
| ``\mathbb{A}`` | ``A_{ij}(\vec{x}, \vec{x}')`` | integral kernel: ``\vert\vec{v}\rangle = \mathbb{A}\vert\vec{u}\rangle`` means ``v_{i}(\vec{x}) = \sum_{j} \int A_{ij}(\vec{x}, \vec{x}')u_{j}(\vec{x}')~\mathrm{d}^{3} x'`` |
| ``\mathbb{I}`` | ``\delta_{ij} \delta^{3}(\vec{x} - \vec{x}')`` | identity operator |
| ``\mathbb{P}_{n}`` | ``\delta_{ij} \delta^{3}(\vec{x} - \vec{x}') \Theta(\vec{x} \in V_{n})`` | projection onto material DOFs of body ``n`` |
| ``\mathbb{A}^{\dagger}`` | ``(\mathbb{A}^{\dagger})_{ij}(\vec{x}, \vec{x}') = A^{\star}_{ji}(\vec{x}', \vec{x})`` | Hermitian adjoint: basis-independent, ``\langle \vec{u}, \mathbb{A}^{\dagger} \vec{v} \rangle = \langle \mathbb{A} \vec{u}, \vec{v} \rangle`` |
| ``\operatorname{sym}(\mathbb{A})`` | ``(\mathbb{A} + \mathbb{A}^{\dagger})/2`` | Hermitian part for square operator |
| ``\operatorname{asym}(\mathbb{A})`` | ``(\mathbb{A} - \mathbb{A}^{\dagger})/(2\mathrm{i})`` | anti-Hermitian part for square operator |
| ``\mathbb{A}^{\top}`` | ``(\mathbb{A}^{\top})_{ij}(\vec{x}, \vec{x}') = A_{ji}(\vec{x}', \vec{x})`` | unconjugated transpose: basis-*dependent* |
| ``\mathbb{A}^{\star}`` | ``(\mathbb{A}^{\star})_{ij}(\vec{x}, \vec{x}') = A^{\star}_{ij}(\vec{x}, \vec{x}')`` | complex conjugate: basis-*dependent* |
| ``\operatorname{Re}(\mathbb{A})`` | ``(\mathbb{A} + \mathbb{A}^{\star})/2`` | real part |
| ``\operatorname{Im}(\mathbb{A})`` | ``(\mathbb{A} - \mathbb{A}^{\star})/(2\mathrm{i})`` | imaginary part |
| ``\mathbb{AB}`` | ``\sum_{l} \int A_{il}(\vec{x}, \vec{x}'') B_{lj}(\vec{x}'', \vec{x}')~\mathrm{d}^{3} x''`` | operator product, *not* commutative
| ``\operatorname{Tr}[\mathbb{A}]`` | ``\sum_{i} \int A_{ii}(\vec{x}, \vec{x})~\mathrm{d}^{3} x`` | trace, basis-independent |

Note also that operators relevant to EM theory will commonly depend on
the frequency ``\omega`` as a parameter; this will generally be
suppressed in the notation for brevity, and context will make clear
whether the dependence is on real frequency ``\omega`` versus
imaginary frequency ``\omega = \mathrm{i}\xi``. Furthermore, the
Maxwell Green's function ``\mathbb{G}`` is defined to satisfy
``[\nabla \times (\nabla \times) - (\omega/c)^{2} (\mathbb{I} +
\mathbb{V})]\mathbb{G} = (\omega/c)^{2} \mathbb{I}`` in the presence
of a (possibly nonlocal, inhomogeneous, anisotropic) susceptibility
``\mathbb{V}``, so that the solution in vacuum (``\mathbb{V} \to 0``)
is ``G^{(0)}_{ij}(\vec{x}, \vec{x}') = [\partial_{i}
\partial_{j} + (\omega/c)^{2}
\delta_{ij}](e^{\mathrm{i}\omega\vert\vec{x} -
\vec{x}'\vert/c}/(4\pi\vert\vec{x} - \vec{x}'\vert))``, and this
reproduces the field due to a static point dipole in the electrostatic
(``c \to \infty``) limit.

## Formulas relevant for FED

AARMBEM.jl assumes that there exist ``N`` molecular bodies labeled
``n`` with susceptibilities ``\mathbb{V}_{n}`` that are
separable/distinct/disjoint from each other (i.e. there are no
inherent nonlocal correlations between bodies); this separability
means the total susceptibility of the system may be written as a sum
``\mathbb{V} = \sum_{n = 1}^{N} \mathbb{V}_{n}``. These interact in
the presence of a macroscopic environment (which, for the
implementation in AARMBEM.jl, may be vacuum or a PEC plane at ``z =
0``) whose Maxwell Green's function is defined as
``\mathbb{G}^{\mathrm{mac}}``; for convenience, the scattering Green's
function will also be defined as ``\mathbb{G}^{\mathrm{sca}} =
\mathbb{G}^{\mathrm{mac}} - \mathbb{G}^{(0)}`` (and this vanishes in
vacuum).

Given these definitions, the most important quantity for EM scattering
in general, and for FED phenomena in turn, is the T-operator, defined
as ``\mathbb{T} = (\mathbb{V}^{-1} - \mathbb{G}^{\mathrm{mac}})^{-1} =
(\mathbb{I} - \mathbb{V}\mathbb{G}^{\mathrm{mac}})^{-1}
\mathbb{V}``. This describes EM scattering to all orders within and
between all bodies involved. It can be defined in block matrix form as
```math
\mathbb{T}^{-1} =
\begin{bmatrix}
\mathbb{V}_{1}^{-1} - \mathbb{G}^{\mathrm{mac}}_{1,1} & -\mathbb{G}^{\mathrm{mac}}_{1,2} & -\mathbb{G}^{\mathrm{mac}}_{1,3} & \ldots \\
-\mathbb{G}^{\mathrm{mac}}_{2,1} & \mathbb{V}_{2}^{-1} - \mathbb{G}^{\mathrm{mac}}_{2,2} & -\mathbb{G}^{\mathrm{mac}}_{2,3} & \ldots \\
-\mathbb{G}^{\mathrm{mac}}_{3,1} & -\mathbb{G}^{\mathrm{mac}}_{3,2} & \mathbb{V}_{3}^{-1} - \mathbb{G}^{\mathrm{mac}}_{3,3} & \ldots \\
\vdots & \vdots & \vdots & \ddots
\end{bmatrix}
```
where the blocks ``\mathbb{G}^{\mathrm{mac}}_{mn} = \mathbb{P}_{m}
\mathbb{G}^{\mathrm{mac}} \mathbb{P}_{n}`` are the projections of
``\mathbb{G}^{\mathrm{mac}}`` onto the DOFs of the molecules ``m`` and
``n``. Also relevant to vdW interactions in particular is the quantity
``\mathbb{T}_{\infty}``, which is the evaluation of ``\mathbb{T}``
when each molecular body is taken to be infinitely far away from every
other molecular and macroscopic body (for which
``\mathbb{G}^{\mathrm{mac}}`` may be replaced by
``\mathbb{G}^{(0)}``).

In terms of this, the vdW interaction free energy at thermal
equilibrium is defined for temperature ``T`` as
```math
\mathcal{F}(T) = k_{\mathrm{B}} T\sum_{l = 0}^{\infty} {}' f(\mathrm{i}\xi_{l})
```
in terms of the integrand ``f(\mathrm{i}\xi) =
\ln(\det[\mathbb{T}_{\infty} \mathbb{T}^{-1}])`` evaluated at the
Matsubara frequencies ``\xi_{l} = 2\pi k_{\mathrm{B}} Tl/\hbar``,
where the prime over the summation implies an additional prefactor of
``1/2`` for the contribution ``l = 0``; at ``T = 0``, this turns into
the integral ``\mathcal{F}(0) = \mathcal{E} = \hbar\int_{0}^{\infty}
f(\mathrm{i}\xi)~\frac{\mathrm{d}\xi}{2\pi}``.

The thermal radiation power may generally be defined as
```math
P = \int_{0}^{\infty} W(\omega)~\frac{\mathrm{d}\omega}{2\pi}
```
in terms of an appropriate function ``W(\omega)``; it is assumed that
each body ``n`` is at a certain temperature ``T_{n}``, and that these
are immersed in an environment of temperature ``T_{\mathrm{env}}``
(which, for the implementation in AARMBEM.jl, is also the temperature
of the PEC plane if present). For total thermal emission from a
specific body ``m``, the function is defined as
```math
W^{(m)} = \sum_{n = 1}^{N} s_{nm} [\Pi(\omega, T_{n}) - \Pi(\omega,
T_{\mathrm{env}})]\Phi^{(m)}_{n}
```
in terms of the sign function ``s_{nm} = 1 - 2\delta_{nm}``, while for
net radiative heat transfer between bodies ``m`` and ``n``, the
function is defined as
```math
W_{m \to n} = [\Pi(\omega, T_{m}) - \Pi(\omega,
T_{n})]\Phi^{(m)}_{n}
```
and both of these are defined in terms of the Planck function
``\Pi(\omega, T) = \hbar\omega/(e^{\hbar\omega/(k_{\mathrm{B}} T)} -
1)`` as well as the spectral function ``\Phi^{(m)}_{n}
(\omega)``. This spectral function is defined as
```math
\Phi^{(m)}_{n} =
-4\operatorname{Tr}[\operatorname{asym}(\mathbb{V}_{m}^{-1\dagger})
\mathbb{P}_{m} \mathbb{T}^{\dagger} \operatorname{asym}(\mathbb{P}_{n}
\mathbb{G}^{\mathrm{mac}}) \mathbb{T} \mathbb{P}_{m}]
```
in terms of the aforementioned macroscopic Green's function and
overall T-operators. Reciprocity means that it satisfies
``\Phi^{(m)}_{n} (\omega) = \Phi^{(n)}_{m} (\omega)``.

These formulas are fully general, and when all bodies are compact
(i.e. have well-defined finite dimensions in all directions
appropriate for the problem at hand), these formulas require no
further analytical manipulation. When all bodies have commensurate
discrete translational symmetries in ``d \in \{1, 2\}`` dimensions, it
is most useful to perform these computations imposing Bloch periodic
boundary conditions in a unit cell of general volume (which in 1
dimension is a length, or in 2 dimensions is an area)
``V_{\mathrm{uc}}``, and where Bloch wavevectors ``\vec{k}`` lie
within the BZ. This means that we use as the vdW free energy integrand
```math
f(\mathrm{i}\xi) = \int_{\mathrm{BZ}} f(\mathrm{i}\xi,
\vec{k})~\frac{V_{\mathrm{uc}}~\mathrm{d}^{d} k}{(2\pi)^{d}}
```
where ``f(\mathrm{i\xi}, \vec{k}) = \ln(\det[\mathbb{T}_{\infty}
\mathbb{T}^{-1}])`` has all relevant operators evaluated at
``(\mathrm{i}\xi, \vec{k})``; this will yield a vdW free energy per
unit cell. Because ``(\mathbb{T}(\mathrm{i}\xi, \vec{k}))^{\top} =
\mathbb{T}(\mathrm{i}\xi, -\vec{k})`` with a similar relation holding
for ``\mathbb{T}_{\infty}``, and the determinant is invariant with
respect to transposition and for a product of operators is invariant
with respect to their ordering, then ``f(\mathrm{i}\xi, \vec{k}) =
f(\mathrm{i}\xi, -\vec{k})``.  Similarly, we use as the thermal
radiation power integrand
```math
\Phi^{(m)}_{n} (\omega) = \int_{\mathrm{BZ}} \Phi^{(m)}_{n} (\omega,
\vec{k})~\frac{V_{\mathrm{uc}}~\mathrm{d}^{d} k}{(2\pi)^{d}}
```
where ``\Phi^{(m)}_{n}(\omega, \vec{k}) =
-4\operatorname{Tr}[\operatorname{asym}(\mathbb{V}_{m}^{-1\dagger})
\mathbb{P}_{m} \mathbb{T}^{\dagger} \operatorname{asym}(\mathbb{P}_{n}
\mathbb{G}^{\mathrm{mac}}) \mathbb{T} \mathbb{P}_{m}]`` has all
relevant operators evaluated at ``(\omega, \vec{k})``; this will yield
a thermal radiation power per unit cell. Because the various operators
are reciprocal, satisfying ``(\mathbb{A}(\omega, \vec{k}))^{\top} =
\mathbb{A}(\omega, -\vec{k})``, then ``\Phi^{(m)}_{n} (\omega,
\vec{k}) = \Phi^{(n)}_{m} (\omega, -\vec{k})``.

The following sections will discuss the details of the RMB framework
in representing molecular polarization response functions to construct
these operators in appropriate bases; these representations can be
generalized to other EM scattering phenomena beyond FED phenomena, but
the focus of AARMBEM.jl is on FED phenomena.

## RMB representation of EM scattering operators

For any of these quantities to be computed, these operators require
representations as matrices. As the molecular bodies are assumed to be
disjoint, the total susceptibility is block-diagonal in the basis of
molecular bodies. Given this, each molecular susceptibility may be
expanded generically in a certain basis. For bodies that have
atom-scale feature sizes but highly delocalized electronic response,
the basis functions should be correspondingly delocalized in real
space. The RMB framework in principle allows for this, but the
AARMBEM.jl code has not implemented this. Instead, the AARMBEM.jl code
implements a set of basis functions appropriate for insulating or
weakly conducting media. This is described as follows.

### Susceptibilities

Appropriate for relatively localized electronic response, each
molecular body has mapped onto itself a set of effective oscillator
DOFs, specifically corresponding to one effective valence electronic
oscillator and one effective screened nuclear oscillator for each
atom. This means the susceptibility of each molecular body may be
written as ``\mathbb{V}_{n} = \sum_{p,i,q,j} \alpha_{pi,qj}
\vert\vec{f}_{pi}\rangle\langle\vec{f}_{qj}\vert``, where the basis
functions ``\vert\vec{f}_{pi}\rangle`` are chosen as appropriate
localized functions, while the expansion coefficients
``\alpha_{pi,qj}`` are the elements of a matrix ``\alpha`` that are
derived from oscillator equations of motion described as follows.

For each atom ``p``, the valence electronic oscillator has an
effective mass ``m_{\mathrm{e}p}``, charge (coupling it to long-range
EM fields) ``q_{\mathrm{e}p}``, damping coefficient
``b_{\mathrm{e}p}``, and isotropic spring constant (coupling it to its
corresponding screened nuclear oscillator) ``k_{\mathrm{e}p}``, while
the screened nuclear oscillator at position ``\vec{r}_{p}`` does not
couple directly to long-range EM fields due to screening by electrons
but does have an effective mass ``m_{\mathrm{I}p}``, damping
coefficient ``b_{\mathrm{I}p}``, and possibly anisotropic spring
constants (coupling it to other nearby screened nuclear oscillators
``q`` along Cartesian directions ``i`` and ``j``)
``K_{\mathrm{I}pi,qj}``. These quantities are derived through a method
briefly summarized as follows. The damping coefficients
``b_{\mathrm{e}p}`` and ``b_{\mathrm{I}p}`` may in principle be
derived ab-initio, but especially in the regime of low loss, it is
sufficient for users to specify them arbitrarily. The effective
nuclear masses ``m_{\mathrm{I}p}`` are taken from elemental data. The
remaining quantities are taken from density functional theory (DFT)
calculations in the following way. The ground state electron density
of the molecular body in isolation is computed using DFT, and the
nuclear coordinates are relaxed to the lowest energy state; the
Hessian of the ground state energy with respect to displacements in
the nuclear spatial coordinates immediately yields the nuclear spring
constants ``K_{\mathrm{I}pi,qj}``, and the nuclear spatial coordinates
``\vec{r}_{p}`` are retained as well. The ground state electron
density then undergoes a Hirshfeld partitioning to give effective
electron densities around each atom, accounting for local chemical
effects relative to the atom in isolation, and this is then combined
with accurate reference data for the static polarizability of the atom
in isolation to give a static polarizability for the atom within the
molecule, which will be called ``\alpha_{\mathrm{e}p0}``. The valence
electronic oscillator has a polarizability
``\alpha_{\mathrm{e}p}(\omega) = \frac{\alpha_{\mathrm{e}p0}}{1 -
(\omega/\omega_{\mathrm{e}p})^{2}}``, with the frequency
``\omega_{\mathrm{e}p}`` taken from accurate reference data in an
Unsoeld approximation. This is mapped to
``\alpha_{\mathrm{e}p}(\omega) =
\frac{q_{\mathrm{e}p}^{2}}{k_{\mathrm{e}p} - \omega^{2}
m_{\mathrm{e}p}}`` with the constraints ``q_{\mathrm{e}p} =
N_{\mathrm{e}p} q_{\mathrm{e}}`` and ``m_{\mathrm{e}p} =
N_{\mathrm{e}p} m_{\mathrm{e}}`` in terms of the known fundamental
electron mass ``m_{\mathrm{e}}`` and charge ``q_{\mathrm{e}}`` to
yield the total number of electrons ``N_{\mathrm{e}p}`` effectively
constituting the valence electronic oscillator. Solving these
equations yields ``q_{\mathrm{e}p} = \alpha_{\mathrm{e}p0}
\omega_{\mathrm{e}p}^{2} m_{\mathrm{e}} / q_{\mathrm{e}}``,
``m_{\mathrm{e}p} = \alpha_{\mathrm{e}p0} \omega_{\mathrm{e}p}^{2}
m_{\mathrm{e}}^{2} / q_{\mathrm{e}}^{2}``, and ``k_{\mathrm{e}p} =
\alpha_{\mathrm{e}p0} \omega_{\mathrm{e}p}^{4} m_{\mathrm{e}}^{2} /
q_{\mathrm{e}}^{2}`` in terms of ``\alpha_{\mathrm{e}p0}`` and
``\omega_{\mathrm{e}p}`` along with fundamental constants.

!!! info

    AARMBEM.jl does *not* perform DFT calculations, nor does it look
    up reference data. All calculations required to compute
    ``\alpha_{\mathrm{e}p0}``, ``\omega_{\mathrm{e}p}``, and
    ``K_{\mathrm{I}pi,qj}`` must be done separately.

    ``K_{\mathrm{I}pi,qj}`` can be obtained with any electronic
    structure code capable of DFT calculations with analytical
    nuclear gradients. The force constants (Hessian) are then obtained
    by numerical finite differencing of the analytical gradients.
    Care must be taken to ensure that the DFT calculations are sufficiiently
    converged to eliminate any numerical noise in the Hessian.

    The Hirshfeld volumes used to scale the reference free-atom polarizability
    can be calculated from the molecular electron density and the respective
    free-atom electron densities, as described in [^TkatchenkoPRL09].

    The reference free-atom polarizabilities and oscillator frequencies can be
    obtained from [^GouldJCTC16].


Once these quantities are computed, they may be arranged into
matrices. The nuclear spring constants ``K_{\mathrm{I}pi,qj}`` can
naturally be arranged into a matrix ``K_{\mathrm{I}}`` due to the
anisotropy and nonlocality of the couplings. The other quantities are
isotropic and local (i.e. they do not involve couplings between
different atoms), so quantities generally of the form ``a_{p}`` can be
assembled as diagonal matrices ``A_{pi,qj} = a_{p} \delta_{pq}
\delta_{ij}``. These are then put together to yield the oscillator
equations of motion
```math
\begin{bmatrix}
K_{\mathrm{e}} - \operatorname{i}\omega B_{\mathrm{e}} - \omega^{2} M_{\mathrm{e}} & -K_{\mathrm{e}} \\
-K_{\mathrm{e}} & K_{\mathrm{e}} + K_{\mathrm{I}} - \operatorname{i}\omega B_{\mathrm{I}} - \omega^{2} B_{\mathrm{I}}
\end{bmatrix}
\begin{bmatrix}
Q_{\mathrm{e}}^{-1} p_{\mathrm{e}} \\
x_{\mathrm{I}}
\end{bmatrix}
=
\begin{bmatrix}
Q_{\mathrm{e}} e_{\mathrm{e}} \\
0
\end{bmatrix}
```
for which the relation ``p_{\mathrm{e}} = \alpha e_{\mathrm{e}}`` is
then used to determine the polarizability matrix ``\alpha`` entering
the susceptibility ``\mathbb{V}_{n}``. Explicitly, the polarizability
matrix is written as
```math
\alpha = Q_{\mathrm{e}} (K_{\mathrm{e}} - \mathrm{i}\omega
B_{\mathrm{e}} - \omega^{2} M_{\mathrm{e}} - K_{\mathrm{e}}
(K_{\mathrm{e}} + K_{\mathrm{I}} - \mathrm{i}\omega B_{\mathrm{I}} -
\omega^{2} M_{\mathrm{I}})^{-1} K_{\mathrm{e}})^{-1} Q_{\mathrm{e}}
```
in terms of the aforementioned matrices.

The effective valence electronic oscillators are taken to be in their
ground states, so with this in mind, the basis functions are chosen to
be the oscillator ground state densities
```math
\vec{f}_{pi}(\vec{x}) = \left(\sqrt{2\pi} \sigma_{p}\right)^{-3}
\exp\left(-\frac{(\vec{x} - \vec{r}_{p})^{2}}{2\sigma_{p}^{2}}\right)
\vec{e}_{i}
```
in each Cartesian direction ``\vec{e}_{i}``. These are defined in
terms of the aforementioned nuclear coordinates ``\vec{r}_{p}``
(arising from DFT) along with the Gaussian widths ``\sigma_{p}``,
which are not arbitrary but are chosen at each frequency to depend on
the polarizability matrix in the following way: the definition
``\sigma_{p}(\omega) = (\alpha_{p}(\omega) / 3)^{1/3} /
(2\sqrt{\pi})`` ensures that the internal self-energy of a
Gaussian-distributed polarization is equal to the energy of its
electrostatic self-interaction, and the definition ``\alpha_{p}
(\omega) = \vert\sum_{q,j} \alpha_{pj,qj} (\omega)\vert/3`` takes the
average over Cartesian polarization directions and transforms the
nonlocal polarizability matrix into a set of effective screened atomic
polarizabilities (with the absolute value ensuring real nonnegative
Gaussian widths even at real ``\omega``).

For compact bodies (i.e. having no dimensions of infinite extent) as
well as for periodic bodies (i.e. having infinite extent with discrete
translational symmetry in 1 or 2 dimensions), the matrices
``Q_{\mathrm{e}}``, ``M_{\mathrm{e}}``, ``B_{\mathrm{e}}``,
``K_{\mathrm{e}}``, ``M_{\mathrm{I}}``, and ``B_{\mathrm{I}}`` are all
diagonal nonnegative matrices (independent of the Bloch wavevector
``\vec{k}`` in the BZ for periodic bodies). For compact bodies,
``K_{\mathrm{I}}`` is a real-symmetric positive-semidefinite matrix,
so ``K_{\mathrm{I}} = K_{\mathrm{I}}^{\top} =
K_{\mathrm{I}}^{\dagger}`` (or, more explicitly, ``K_{\mathrm{I}pi,qj}
= K_{\mathrm{I}pi,qj}^{\star} = K_{\mathrm{I}qj,pi}``). This yields a
reciprocal (but *not Hermitian, except at imaginary frequency* where
it is real-symmetric positive-definite) polarizability matrix
satisfying ``\alpha = \alpha^{\top}`` (or, more explicitly,
``\alpha_{pi,qj} = \alpha_{qj,pi}``). For periodic bodies, evaluation
of ``K_{\mathrm{I}\vec{k}}`` must be done for each Bloch wavevector
``\vec{k}`` in the BZ, and ``K_{\mathrm{I}\vec{k}}`` is a Hermitian
positive-semidefinite matrix, so ``K_{\mathrm{I}\vec{k}} =
(K_{\mathrm{I}\vec{k}})^{\dagger}``
(i.e. ``(K_{\mathrm{I}\vec{k}})^{\dagger}_{pi,qj} =
K^{\star}_{\mathrm{I}\vec{k}qj,pi}``), while reciprocity means that
``(K_{\mathrm{I}\vec{k}})^{\top} = K_{\mathrm{I}(-\vec{k})}``
(i.e. ``(K_{\mathrm{I}\vec{k}})^{\top}_{pi,qj} =
K_{\mathrm{I}(-\vec{k})qj,pi}``). This yields the reciprocity relation
for the polarizability matrix of a periodic body, namely that
``(\alpha_{\vec{k}})^{\top} = \alpha_{-\vec{k}}``
(i.e. ``(\alpha_{\vec{k}})^{\top}_{pi,qj} =
\alpha_{(-\vec{k})qj,pi}``), while specifically at imaginary frequency
``\omega = \mathrm{i}\xi``, ``\alpha_{\vec{k}}`` is a Hermitian
positive-definite matrix, so ``\alpha_{\vec{k}} =
\alpha_{\vec{k}}^{\dagger}``.

### Green's functions

Once the polarizability matrix is constructed, all that is left to
represent the operators necessary for computing FED and other EM
quantities is to represent the Green's functions in matrix form. In
particular, the matrix elements ``\langle \vec{f}_{pi},
\mathbb{G}^{\mathrm{mac}} \vec{f}_{qj} \rangle`` are required. For a
general macroscopic environment, this almost always requires costly
six-dimensional cubature, and that is in the *easier* case of
semianalytical expression of ``\mathbb{G}^{\mathrm{mac}}``. For this
reason, the only macroscopic environments currently supported in
AARMBEM.jl are either vacuum or a single PEC plane at ``z = 0``,
because the matrix elements ``\langle \vec{f}_{pi}, \mathbb{G}^{(0)}
\vec{f}_{qj} \rangle`` can be expressed analytically, and image theory
allows using those expressions to determine ``\langle \vec{f}_{pi},
\mathbb{G}^{\mathrm{mac}} \vec{f}_{qj} \rangle`` in the presence of a
PEC plane at ``z = 0``.

#### Compact media

For compact media, the matrix elements are computed directly in real
space. In particular, they are written as
```math
G^{(0)}_{pi,qj} = \langle \vec{f}_{pi}, \mathbb{G}^{(0)} \vec{f}_{qj}
\rangle = \left[\frac{\partial^{2}}{\partial r_{pi} \partial r_{pj}} +
(\omega/c)^{2} \delta_{ij}\right] \frac{\exp(-q^{2} /
4)}{8\pi\vert\vec{r}_{p} - \vec{r}_{q}\vert} \left[e^{\mathrm{i}\rho
q} \operatorname{erfc} \left(-\frac{\mathrm{iq}}{2} - \rho\right) -
e^{-\mathrm{i}\rho q} \operatorname{erfc} \left(-\frac{\mathrm{iq}}{2}
+ \rho\right)\right]
```
in terms of the complementary error function ``\operatorname{erfc}(u)
= \frac{2}{\sqrt{\pi}} \int_{u}^{\infty} e^{-v^{2}}~\mathrm{d}v``, and
in terms of the definitions of the dimensionless variables ``q \equiv
(\omega/c)\sqrt{2(\sigma_{p}^{2} + \sigma_{q}^{2})}`` and ``\rho
\equiv \vert\vec{r}_{p} - \vec{r}_{q}\vert/\sqrt{2(\sigma_{p}^{2} +
\sigma_{q}^{2})}``. The diagonal elements ``G^{(0)}_{pi,pj}``
corresponding to the self-interaction of an oscillator in vacuum are
taken to vanish when contribution to the T-operators for vdW
interaction free energies or thermal radiation powers, but not in the
contribution of the term ``\operatorname{asym}(\mathbb{P}_{n}
\mathbb{G}^{\mathrm{mac}})`` to the thermal radiation power; the
diagonal elements ``G^{\mathrm{sca}}_{pj,pj}`` in the presence of a
PEC plane are never taken to vanish, as the corresponding vacuum
interaction is with a physically separated image dipole, not with the
original dipole in itself. All of these quantities may be evaluated at
real ``\omega`` or at ``\omega = \mathrm{i}\xi`` for real nonnegative
``\xi``. At real frequency ``\omega``, the matrices ``G^{(0)}`` and
``G^{\mathrm{mac}}`` are reciprocal, meaning ``G^{(0)} = G^{(0)\top}``
and ``G^{\mathrm{mac}} = G^{\mathrm{mac}\top}``; at imaginary
frequency ``\omega = \mathrm{i}\xi``, these matrices are real-valued
as well.

#### Periodic media

For periodic media, the matrix elements ``G^{(0)}_{\vec{k}pi,qj} =
\sum_{\vec{R}} e^{-\mathrm{i} \vec{k} \cdot \vec{R}} \langle \vec{f}_{p + \vec{R}, i}, \mathbb{G}^{(0)} \vec{f}_{qj} \rangle`` are needed in
Bloch space, with Bloch wavevector ``\vec{k}`` restricted to lie in
the BZ; here, ``\vert\vec{f}_{p+\vec{R},i}\rangle`` has position space
representation ``\vec{f}_{pi}(\vec{x} - \vec{R})``, where ``\vec{R}``
refers to a lattice vector. In any periodic dimensionality, it is
useful to write ``G^{(0)}_{\vec{k}pi,qj} =
G^{(0)\mathrm{SR}}_{\vec{k}pi,qj} +
G^{(0)\mathrm{LR}}_{\vec{k}pi,qj}`` and make use of the Ewald
summation procedure to efficiently contribute each of these two terms
contributing to the sum. Further, in any periodic dimensionality,
reciprocity means that ``(G^{(0)}_{\vec{k}})^{\top} =
G^{(0)}_{-\vec{k}}`` at any frequency ``\omega``, while at imaginary
frequency ``\omega = \mathrm{i}\xi``, ``(G^{(0)}_{\vec{k}})^{\dagger}
= G^{(0)}_{\vec{k}}``; these relations also hold for
``G^{(0)\mathrm{SR}}_{\vec{k}}`` and ``G^{(0)\mathrm{LR}}_{\vec{k}}``,
and for ``G^{\mathrm{mac}}_{\vec{k}}`` in the presence of a PEC plane
at ``z = 0``.

For either 1 or 2 periodic dimensions, given lattice vectors
``\vec{R}``, the first contribution may be written at imaginary
frequency ``\omega = \mathrm{i}\xi`` as
```math
G^{(0)\mathrm{SR}}_{\vec{k}pi,qj} = \left[\frac{\partial^{2}}{\partial r_{pi} \partial r_{pj}} - (\xi/c)^{2} \delta_{ij}\right] \times \sum_{\vec{R}} \frac{e^{(\sigma_{p}^{2} + \sigma_{q}^{2})\xi^{2} / (2c^{2}) - \mathrm{i}\vec{k} \cdot \vec{R}}}{8\pi\vert\vec{r}_{p} + \vec{R} - \vec{r}_{q}\vert} \times \\
\left\{e^{-\rho\theta} \left[\operatorname{erfc}\left(\nu\rho - \frac{\theta}{2\nu}\right) - \operatorname{erfc}\left(\rho - \frac{\theta}{2}\right)\right] + e^{\rho\theta} \left[\operatorname{erfc}\left(\nu\rho + \frac{\theta}{2\nu}\right) - \operatorname{erfc}\left(\rho + \frac{\theta}{2}\right)\right]\right\}
```
having defined the dimensionless variables ``\rho \equiv
\vert\vec{r}_{p} + \vec{R} - \vec{r}_{q}\vert/\sqrt{2(\sigma_{p}^{2} +
\sigma_{q}^{2})}`` and ``\theta \equiv \sqrt{2(\sigma_{p}^{2} +
\sigma_{q}^{2})}\xi/c``, and where ``\nu \in [0, 1]`` is a
dimensionless parameter that is ultimately chosen by the user for
convergence. The expression at real frequency ``\omega`` can be
obtained be consistently substituting ``\xi = -\mathrm{i}\omega``.

For 1 periodic dimension, it shall be assumed that all lattice vectors
can be written as ``\vec{R} = n\vec{a}`` for all integers ``n`` given
a primitive lattice vector ``\vec{a} = a\vec{e}_{\vec{a}}`` (such that
``a = \vert\vec{a}\vert`` and ``\vec{e}_{\vec{a}} = \vec{a}/a``), and
the equation ``\vec{r}_{p} - \vec{r}_{q} = \Delta r_{\parallel}
\vec{e}_{\vec{a}} + \Delta \vec{r}_{\perp}`` shall define ``\Delta
r_{\parallel}`` and ``\Delta \vec{r}_{\perp}``. Furthermore, the
primitive reciprocal lattice vector is ``\vec{b} =
2\pi\vec{e}_{\vec{a}}/a``, reciprocal lattice vectors will satisfy
``\vec{g} = n\vec{b}`` for all integers ``n``, and ``\vec{k}`` in the
BZ will lie along (or opposite to) ``\vec{e}_{\vec{a}}``. At imaginary
frequency ``\omega = \mathrm{i}\xi``, this leads to the expression
```math
G^{(0)\mathrm{LR}}_{\vec{k}pi,qj} = \frac{|\vec{b}|}{8\pi^{2}} \left[\frac{\partial^{2}}{\partial r_{pi} \partial r_{pj}} - (\xi/c)^{2} \delta_{ij}\right] \times \\
\sum_{\vec{g}} \left[\exp(\theta^{2} / 4 + \mathrm{i}(\vec{k} + \vec{g}) \cdot (\vec{r}_{p} - \vec{r}_{q})) \sum_{s = 0}^{\infty} \frac{(-1)^{s}}{s!} (\nu\rho_{\perp})^{2s} E_{s + 1} (\eta^{2} / (4\nu^{2}))\right]
```
written again in terms of ``\theta \equiv \sqrt{2(\sigma_{p}^{2} +
\sigma_{q}^{2})}\xi/c``, and where ``\nu \in [0, 1]`` is a
dimensionless parameter that is ultimately chosen by the user for
convergence, and also in terms of new dimensionless quantities
``\rho_{\perp} =
\vert\Delta\vec{r}_{\perp}\vert/\sqrt{2(\sigma_{p}^{2} +
\sigma_{q}^{2})}`` and ``\eta^{2} = \theta^{2} + 2(\sigma_{p}^{2} +
\sigma_{q}^{2})\vert\vec{k} + \vec{g}\vert^{2}``. Also used is the
definition of the exponential integral function ``E_{s + 1}(x) =
s^{-1} (e^{-x} - xE_{s}(x))`` given ``E_{1}(x) = \int_{x}^{\infty}
t^{-1} e^{-t}~\mathrm{d}t``. Once again, the expression at real
frequency ``\omega`` can be obtained be consistently substituting
``\xi = -\mathrm{i}\omega``.

For 2 periodic dimensions, it shall be assumed that all lattice
vectors can be written as ``\vec{R} = n_{1} \vec{a}_{1} + n_{2}
\vec{a}_{2}`` for all integer pairs ``(n_{1}, n_{2})``, even if the
primitive lattice vectors ``\vec{a}_{1}`` and ``\vec{a}_{2}`` are
neither orthogonal nor unit-normalized. Furthermore, if the primitive
lattice vectors are orthonormalized to yield a basis ``\vec{e}_{1}``
and ``\vec{e}_{2}`` defining the plane of periodicity and normal
vector ``\vec{e}_{\perp}``, then this allows for writing reciprocal
lattice vectors ``\vec{g} = n_{1} \vec{b}_{1} + n_{2} \vec{b}_{2}``
for all integer pairs ``(n_{1}, n_{2})`` in terms of the primitive
reciprocal lattice vectors ``\vec{b}_{1} = \frac{2\pi \vec{a}_{2}
\times \vec{e}_{\perp}}{\vec{a}_{1} \cdot (\vec{a}_{2} \times
\vec{e}_{\perp})}`` and ``\vec{b}_{2} = \frac{2\pi \vec{e}_{\perp}
\times \vec{a}_{1}}{\vec{a}_{2} \cdot (\vec{e}_{\perp} \times
\vec{a}_{1})}``, and defining ``\vec{r}_{p} - \vec{r}_{q} = \Delta
r_{1} \vec{e}_{1} + \Delta r_{2} \vec{e}_{2} + \Delta r_{\perp}
\vec{e}_{\perp}``. At imaginary frequency ``\omega = \mathrm{i}\xi``,
this leads to the expression
```math
G^{(0)\mathrm{LR}}_{\vec{k}pi,qj} = \frac{\vert\vec{b}_{1} \times \vec{b}_{2}\vert}{16\pi^{2}} \left[\frac{\partial^{2}}{\partial r_{pi} \partial r_{pj}} - (\xi/c)^{2} \delta_{ij}\right] \times \\
\sum_{\vec{g}} \eta^{-1} \sqrt{2(\sigma_{p}^{2} + \sigma_{q}^{2})} \exp(\theta^{2} / 4 + \mathrm{i}(\vec{k} + \vec{g}) \cdot (\vec{r}_{p} - \vec{r}_{q})) \left[e^{-\eta\rho_{\perp}} \operatorname{erfc}\left(\frac{\eta}{2\nu} - \nu\rho_{\perp}\right) + e^{\eta\rho_{\perp}} \operatorname{erfc}\left(\frac{\eta}{2\nu} + \nu\rho_{\perp}\right)\right]
```
written again in terms of
``\theta \equiv \sqrt{2(\sigma_{p}^{2} + \sigma_{q}^{2})}\xi/c`` and
``\eta^{2} = \theta^{2} + 2(\sigma_{p}^{2} + \sigma_{q}^{2})\vert\vec{k} + \vec{g}\vert^{2}``,
and where ``\nu \in [0, 1]`` is a dimensionless parameter that is
ultimately chosen by the user for convergence, and also in terms of new
dimensionless quantity
``\rho_{\perp} = \Delta r_{\perp}/\sqrt{2(\sigma_{p}^{2} + \sigma_{q}^{2})}``
(which now may be negative). Once again, the
expression at real frequency ``\omega`` can be obtained be
consistently substituting ``\xi = -\mathrm{i}\omega``.

!!! note

    In the above expressions, ``\nu \in [0, 1]`` must be consistently
    chosen and applied to computation of both ``G^{(0)\mathrm{SR}}``
    and ``G^{(0)\mathrm{LR}}``. Additionally, the implementation in
    AARMBEM.jl asks not for ``\nu`` but for a related quantity.

!!! warning

    Convergence is quite sensitive to the choice of ``\nu``, which
    will depend heavily on the geometric as well as material
    properties of the molecular body.

[^VenkataramARXIV2020]: Prashanth S. Venkataram, Jan Hermann, Alexandre Tkatchenko, and Alejandro W. Rodriguez. "Fluctuational Electrodynamics in Atomic and Macroscopic Systems: van der Waals Interactions and Radiative Heat Transfer". [arXiv:2005.04083](https://arxiv.org/abs/2005.04083)
[^TkatchenkoPRL09]: Alexandre Tkatchenko, Matthias Scheffler. "Accurate Molecular Van Der Waals Interactions from Ground-State Electron Density and Free-Atom Reference Data". Phys. Rev. Lett. 102, 073005 (2009) [DOI:10.1103/PhysRevLett.102.073005](http://dx.doi.org/10.1103/PhysRevLett.102.073005)
[^GouldJCTC16]: Tim Gould, Tomáš Bučko. "C6 Coefficients and Dipole Polarizabilities for All Atoms and Many Ions in Rows 1–6 of the Periodic Table". J. Chem. Theory Comput. 12, 3603 (2016). [DOI:10.1021/acs.jctc.6b00361](http://dx.doi.org/10.1021/acs.jctc.6b00361)
