# Output quantities in AARMBEM.jl

AARMBEM.jl only outputs dimensionless integrands. It does not perform
any integrations over frequency (or over Bloch wavevector, if
present); this should be done by the user by appropriate sampling and
interpolation over frequency (and Bloch wavevector, if appropriate).

For compact molecular bodies, AARMBEM.jl only requires the
frequency. For vdW interaction free energies, it outputs the integrand
``f(\mathrm{i}\xi)``. For thermal radiation powers, it outputs
one-fourth of the spectrum, namely ``\frac{1}{4} \Phi^{(m)}_{n}
(\omega)``, without any factors of the Planck function.

For molecular bodies with commensurate discrete translational
symmetries in ``d \in \{1, 2\}`` dimensions defining a unit cell of
generalized volume ``V_{\mathrm{uc}}`` and a BZ, AARMBEM.jl requires
the frequency as well as the Bloch wavevector ``\vec{k}``. For vdW
interaction free energies, it outputs the integrand ``f(\mathrm{i}\xi,
\vec{k})``. For thermal radiation powers, it outputs one-fourth of the
spectrum, namely ``\frac{1}{4} \Phi^{(m)}_{n} (\omega, \vec{k})``,
without any factors of the Planck function.

AARMBEM.jl does not natively support computation of vdW forces or
torques. These can be computed by appropriate setup of geometrical
transformations, followed by numerical interpolation and
differentiation with respect to the appropriate distance, angle, or
other parameter.