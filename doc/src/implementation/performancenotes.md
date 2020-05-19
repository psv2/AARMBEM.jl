# Performance notes and required packages

AARMBEM.jl has been written in Julia. It has been tested in Julia
1.0. It requires the modules `LinearAlgebra`, `DelimitedFiles`,
`SharedArrays`, `SpecialFunctions`, `Distributed`, and `Printf`.

AARMBEM.jl has been written with parallelism in mind. In particular,
the functions for assembling the Green's function matrices make use of
parallel for-loops. However, while the programs in AARMBEM.jl may be
invoked using parallel processing via a command like `julia -p N`
(followed by the command as well as any arguments) for a number of
processors `N` larger than 1, testing thus far has showed this to be
somewhat unreliable on certain setups. Therefore, it is recommended
that `N` be set to 1, and that parallelism be done instead over
frequencies (and Bloch wavevectors, if appropriate), as that process
can be done in an embarrassingly parallel manner (i.e. with each set
of frequencies (and Bloch wavevectors, if appropriate) constituting a
single job). This is particularly suitable for use in a script,
especially one that interfaces with a job scheduler like SLURM.

Even with embarrassing parallelism (using one frequency, and one Bloch
wavevector, if appropriate, per job, and one processor per job),
AARMBEM.jl has been seen to occasionally fail to produce output in
certain random instances. These generally have not been predictable,
so users who wish to run a calculation over a large number of
frequencies (and Bloch wavevectors, if appropriate) for a large number
of geometric configurations of a given set of molecular and
macroscopic bodies may need to rerun calculations for certain
combinations of frequency (and Bloch wavevector, if appropriate) and
geometric configuration after the fact.

As discussed in other pages, the code for molecular bodies of infinite
extent obeying Bloch periodic boundary conditions is sometimes less
numerically stable. As an example, it has been observed that
calculations of the vdW interaction free energy integrand of a
graphene sheet at a distance ``z`` above and parallel to a PEC plane
(at ``z = 0``) for a few combinations ``(\mathrm{i}\xi, \vec{k}, z)``
(far fewer than the total number sampled) has the wrong sign or
magnitude. Such occurrences are rare and isolated incidents and can be
mitigated through interpolation as needed, though there may be other
systems where none of the data is reliable or trustworthy; users must
use their best judgment to determine how to proceed.