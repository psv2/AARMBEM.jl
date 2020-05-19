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

In the future, AARMBEM.jl may be properly packaged as a module that
may be installed and used. For now, the simplest way to install it is
to download and locally save the full directory of code from GitHub,
and running the code from other directories pointing to this one.