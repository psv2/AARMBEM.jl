# Ewald summation in AARMBEM.jl

!!! warning

    The implementation of Ewald summation for the matrix elements
    ``G^{(0)}_{\vec{k}pi,qj}`` has several issues of which users
    should be aware.
    
    1. The code has been tested somewhat more for 2 periodic dimensions than for 1 periodic dimension, but both, particularly for 1 periodic dimension, could use more testing.
    
    2. The sum over the order ``s`` of the exponential integral function in computing ``G^{(0)\mathrm{LR}}_{\vec{k}pi,qj}`` in 1 periodic dimension may lead to slow performance depending on the system.

    3. In both 1 and 2 periodic dimensions, the sum over real lattice vectors ``\vec{R}`` when computing ``G^{(0)\mathrm{SR}}_{\vec{k}pi,qj}`` and reciprocal lattice vectors ``\vec{g}`` when computing ``G^{(0)\mathrm{LR}}_{\vec{k}pi,qj}`` must be truncated at some point, and the value of the parameter ``\nu \in [0, 1]`` must be chosen appropriately. In AARMBEM.jl, **these are hardcoded** into the respective functions, in 1 dimension for convergence for a carbyne wire, in 2 dimensions for convergence for a graphene sheet. Consideration of other materials may likely require changing these hardcoded values.

    4. The hardcoded convergence parameter choices in 2 periodic dimensions have been optimized for the RMB model of an infinite sheet of graphene. However, it is not guaranteed that other choices will actually work. For instance, the RMB model of an infinite sheet of single-layer hexagonal BN has so far failed to produce numerically well-behaved results for the ``G^{(0)}_{\mathrm{k}}`` and ``G^{\mathrm{mac}}_{\vec{k}}`` matrices at ``\omega = \mathrm{i}\xi``.

    Given all of these issues, which stem from more general practical
    issues with Ewald summation of ``\mathbb{G}^{(0)}`` even with
    other basis functions, users are advised, when possible, to
    consider avoiding treatment of infinite periodic molecular bodies
    and instead consider using sufficiently large supercells of finite
    size (which would still count as compact bodies). By image theory,
    all of the above issues and warnings also apply to the use of
    Ewald summation to compute ``G^{\mathrm{mac}}_{\vec{k}pi,qj}``
    above a PEC plane at ``z = 0``.