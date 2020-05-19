# Command-line arguments

The two main programs in AARMBEM.jl are AARMBEM_heat.jl, computing
thermal radiation spectra, and AARMBEM_vdW.jl, computing vdW
interaction free energy integrands. Their command-line arguments are
very similar, so this page explains how to run each.

For each, the basic syntax of running the program is generally
```
julia [-p N] PROGRAM arg1=val1 arg2=val2 [...]
```
where the term `[-p N]` refers to the optional argument of running
with a specific number of processors `N`, `PROGRAM` is either
`AARMBEM_heat.jl` or `AARMBEM_vdW.jl`, and `arg1`, `arg2`, and so on
are command-line argument keywords, with `val1`, `val2`, and so on
being the actual corresponding argument values. The ordering of
arguments does not matter, but there should be exactly one space
between the end of one argument value and the beginning of the next
argument keyword, and there should be no spaces between the `=` and
the end of the argument keyword or the beginning of the argument
value.

## Arguments with values

The arguments of the form `arg=val` are as follows.

| `arg` | `val` | Mandatory? | Default behavior (in absence) |
|:------|:------|:-----------|:------------------------------|
| `mollistfilename` | Name of molecular configuration file | Yes | Error |
| `outfilename` | Name of output file for calculations | Yes | Error |
| `translistfilename` | Name of molecular transformation file | No, but strongly recommended | All centers of mass at origin, all orientations preserved |
| `periodicfilename` | Name of Bloch periodicity file | No, but strongly recommended for Bloch periodic systems | Irrelevant for compact geometries, effective supercell for periodic geometries |
| `Genv` | `PEC` if above PEC plane, `VAC` if in vacuum | No | Vacuum (even if the keyword and `val` are present but `val` is anything other than exactly `PEC`) |
| `FloatType` | Julia floating-point data type | No | `Float64` |
| `freqlistfilename` | Name of frequency list file | One of `freqlistfilename` or `freq` is required | Check for `freq` argument keyword |
| `freq` | Floating-point value (in units of radians per second) for single frequency | One of `freqlistfilename` or `freq` is required | Error |
| `klistfilename` | Name of Bloch wavevector list file | Only for Bloch periodicity, at least one of `klistfilename`, `kx`, `ky`, or `kz` is required (`klistfilename` supersedes others) | Check for `kx`, `ky`, or `kz` |
| `kx` | Floating-point value (in units of radians per meter) for single ``k_{x}`` | Only for Bloch periodicity, at least one of `klistfilename`, `kx`, `ky`, or `kz` is required (`klistfilename` supersedes others) | Check for `ky` or `kz` |
| `ky` | Floating-point value (in units of radians per meter) for single ``k_{y}`` | Only for Bloch periodicity, at least one of `klistfilename`, `kx`, `ky`, or `kz` is required (`klistfilename` supersedes others) | Check for `kz` |
| `kz` | Floating-point value (in units of radians per meter) for single ``k_{z}`` | Only for Bloch periodicity, at least one of `klistfilename`, `kx`, `ky`, or `kz` is required (`klistfilename` supersedes others) | Error |

## Arguments without values

There are also arguments that are simply specified as `arg`, as
follows. These arguments are all optional.

Using the argument `nophonons` replaces all polarizability matrices
``\alpha_{n}``, and corresponding atomic polarizabilities entering the
Gaussian basis functions, with purely electronic polarizabilities
``\alpha_{\mathrm{e}n} = Q_{\mathrm{e}n} (K_{\mathrm{e}n} -
\mathrm{i}\omega B_{\mathrm{e}n} - \omega^{2} M_{\mathrm{e}n})^{-1}
Q_{\mathrm{e}n}``. In its absence, the polarizability matrices
``\alpha_{n}`` are computed as usual accounting for the effects of
screened nuclear oscillators. Additionally, the value for the argument
`outfilename` is altered: for example, the value `myoutfile.out` will
be changed to `nophonons_outfilename.out`.

Using the argument `nonretarded` evaluates all Green's function
matrices at ``\omega = 0``, effectively reducing all long-range EM
interactions to Coulomb interactions (whether in vacuum or in the
presence of a PEC plane). In its absence, the Green's function
matrices are evaluated at each ``\omega`` (whether real or imaginary)
as usual. Additionally, the value for the argument `outfilename` is
altered: for example, the value `myoutfile.out` will be changed to
`nonretarded_outfilename.out`, and if the argument `nophonons` is also
present, the value will be changed to
`nonretarded_nophonons_outfilename.out`.