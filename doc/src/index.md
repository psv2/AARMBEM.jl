# AARMBEM.jl

AARMBEM.jl is a code to compute electromagnetic (EM) interactions,
particularly fluctuational EM interactions, at the nanoscale, via the
retarded many-body (RMB) framework. In particular, the RMB framework
computes fluctuational EM interactions among material bodies that may
scale from single atoms or small compact molecules to infinitely long
atomic-scale wires or sheets, in vacuum or in the presence of
macroscopic bodies. AARMBEM.jl is a computational implementation of
the RMB framework.

Currently, two codes are available. One computes van der Waals (vdW)
interaction energies in such systems, while the other computes thermal
radiation powers in such systems. Thus, the focus is on fluctuational
phenomena. However, the theory and the API are general enough that
others may extend this code to perform computations of deterministic
EM phenomena, like absorbed or scattered powers from a specified
incident field/polarization source, local densities of states, et
cetera.

## Name

The name "AARMBEM" is an acronym, expanded to "Ab-initio Atomistic
Retarded Many-Body Electromagnetics at the Mesoscale". It may be
pronounced as "aarambham", identical to a word common in Indian
languages meaning "beginning" or "ab-initio".

## Installation

AARMBEM.jl is a proper Julia package and can be installed as such.
For a quick start, follow these steps:

1. Create a new Julia environment with

```
$ mkdir Project && cd Project
$ julia   # follow by pressing `]`
pkg> activate .
```

2. Install AARMBEM.jl with

```
pkg> add https://github.com/psv2/AARMBEM.jl.git
```

3. Exit package manager with `CTRL+C` and import the `AARMBEM` package:

```
julia> import AARMBEM
```

The command-line programs that use the `AARMBEM` package, `AARMBEM_heat.jl` and
`AARMBEM_vdW.jl`, can be downloaded from
[Github](https://github.com/psv2/AARMBEM.jl) in the subdirectory `scripts`.
