# Output file formats

The code AARMBEM_vdW.jl outputs, to the chosen output file, the vdW
interaction free energy integrand. Each line has the following format:
```
translabel freq f
```
where for compact molecular bodies, `translabel` is one of the
transformation label strings in the molecular transformation files,
`freq` is frequency specified as a real floating-point number
(corresponding to ``\xi``) in scientific notation, and `f` is the
dimensionless integrand ``f(\mathrm{i}\xi)`` for that combination of
geometrical transformation and imaginary frequency. For molecular
bodies with Bloch periodicity in at least one dimension, the output
format is similar, with each line having the following format:
```
translabel freq kx ky kz f
```
where the only differences compared to above are the presence of the
real numbers `kx`, `ky`, and `kz` given in scientific notation
corresponding to the Bloch wavevector components
``(k_{x}, k_{y}, k_{z})``, and that `f` now refers to
``f(\mathrm{i}\xi, \vec{k})``. The integrand value `f` should always
be negative, and any deviations of this would be due to numerical
error, hopefully restricted to places where the integrand is
essentially zero anyway.

The code AARMBEM_heat.jl outputs, to the chosen output file, the
thermal radiation spectrum. For `N` compact molecular bodies, each
line has the following format:
```
translabel freq f[1, 1] f[1, 2] [...] f[1, N] f[2, 1] f[2, 2] [...] f[2, N] [...] f[N, 1] f[N, 2] [...] f[N, N]
```
where now each entry `f[m, n]` is the real floating-point number equal
to the dimensionless quantity ``\frac{1}{4} \Phi^{(m)}_{n} (\omega)``
for that combination of geometrical transformation and frequency, and
where ``\Phi^{(m)}_{n} (\omega) = \Phi^{(n)}_{m} (\omega)``. The
entries `f[n, n]` should be nonpositive, while the entries `f[m, n]`
for `m != n` should be nonnegative, and any deviations of this would
be due to numerical error, hopefully restricted to places where the
integrand is essentially zero anyway.. Similarly, for `N` molecular
bodies with Bloch periodic boundary conditions, each line has the
following format:
```
translabel freq kx ky kz f[1, 1] f[1, 2] [...] f[1, N] f[2, 1] f[2, 2] [...] f[2, N] [...] f[N, 1] f[N, 2] [...] f[N, N]
```
where now each entry `f[m, n]` is the real floating-point number equal
to the dimensionless quantity ``\frac{1}{4} \Phi^{(m)}_{n} (\omega,
\vec{k})`` for that combination of geometrical transformation and
frequency, and where ``\Phi^{(m)}_{n} (\omega, \vec{k}) =
\Phi^{(n)}_{m} (\omega, -\vec{k})``. The entries `f[n, n]` should be
nonpositive, while the entries `f[m, n]` for `m != n` should be
nonnegative, and any deviations of this would be due to numerical
error, hopefully restricted to places where the integrand is
essentially zero anyway..