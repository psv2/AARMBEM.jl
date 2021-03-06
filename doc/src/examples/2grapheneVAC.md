# Two identical parallel graphene sheets (vacuum)

This page links to examples of input files needed to compute vdW
interaction free energies and thermal radiation powers in a system of
Bloch periodic molecular bodies. The particular system under
consideration is of two parallel graphene sheets in vacuum at a
separation ``z`` from each other (with no horizontal offset); each
sheet is identical, and a ``3.9~\mathrm{nm} \times 3.4~\mathrm{nm}``
rectangular unit cell in the ``xy``-plane (i.e. not the hexagonal
primitive unit cell) is used.

!!! note

    The Ewald summation parameters have been hardcoded to give good
    convergence for graphene, as in this example.

## Molecular configuration files

The input data files are as follows.

``B_{\mathrm{e}}`` for graphene: [graphene3p9nmby3p4nmunitcellBe.txt](graphene3p9nmby3p4nmunitcellBe.txt)

``B_{\mathrm{I}}`` for graphene: [graphene3p9nmby3p4nmunitcellBI.txt](graphene3p9nmby3p4nmunitcellBI.txt)

``K_{\mathrm{e}}`` for graphene: [graphene3p9nmby3p4nmunitcellKe.txt](graphene3p9nmby3p4nmunitcellKe.txt)

``K_{\mathrm{I}}`` for graphene: [graphene3p9nmby3p4nmunitcellKI.txt](graphene3p9nmby3p4nmunitcellKI.txt)

``M_{\mathrm{e}}`` for graphene: [graphene3p9nmby3p4nmunitcellMe.txt](graphene3p9nmby3p4nmunitcellMe.txt)

``M_{\mathrm{I}}`` for graphene: [graphene3p9nmby3p4nmunitcellMI.txt](graphene3p9nmby3p4nmunitcellMI.txt)

``Q_{\mathrm{e}}`` for graphene: [graphene3p9nmby3p4nmunitcellQe.txt](graphene3p9nmby3p4nmunitcellQe.txt)

Atomic coordinates ``\{\vec{r}_{p}\}`` for graphene: [graphene3p9nmby3p4nmunitcellxyz.txt](graphene3p9nmby3p4nmunitcellxyz.txt)


Overall configuration file: [2graphene3p9nmby3p4nmVACconfig_periodic.txt](2graphene3p9nmby3p4nmVACconfig_periodic.txt)

## Transformation file

Transformation file: [2graphene3p9nmby3p4nmVACtrans_periodic.txt](2graphene3p9nmby3p4nmVACtrans_periodic.txt)

## Bloch periodicity file

Bloch periodicity file: [1graphene3p9nmby3p4nmlatticevecs.txt](1graphene3p9nmby3p4nmlatticevecs.txt)

## Frequency and Bloch wavevector files

Frequency file: [freq31.txt](freq31.txt)
Bloch wavevector file: [graphene3p9nmby3p4nmunitcellklist1.txt](graphene3p9nmby3p4nmunitcellklist1.txt)

(Many more frequencies must be used to replicate the output for vdW
interaction free energies and thermal radiation powers: these can in
principle be extracted from the corresponding output files themselves,
but this is hampered in practice by the truncation of the output files
in this example for space reasons. For clarity,
in the example output files used, for vdW interactions, the input
frequencies ``\omega = \mathrm{i}|w|`` used were
``w_{n} = 10^{10 + (n - 1)/10}~\mathrm{rad/s}`` for
``n \in \{1, 2, \ldots, 81\}``, while for thermal radiation, the input
frequencies ``\omega = |w|`` used
were ``w_{n} = n \times 10^{12}~\mathrm{rad/s}`` for
``n \in \{1, 2, \ldots, 500\}``. Similarly, many more Bloch wavevectors
must be used. These have been sampled in a less regular manner than
the frequencies but are sampled in the same way for both vdW
interactions and thermal radiation, though the rotational and mirror
symmetries of the rectangular unit cell mean that only
``k_{x} \geq 0`` and ``k_{y} \geq 0`` need to be sampled.)

## Output files

vdW interaction free energies: [2graphene3p9nmby3p4nmVACvdW_periodic.out](2graphene3p9nmby3p4nmVACvdW_periodic.out)

Thermal radiation powers: [2graphene3p9nmby3p4nmVACheat_periodic.out](2graphene3p9nmby3p4nmVACheat_periodic.out)

!!! note

    The ordering of output lines may vary from one run to another.

!!! note

    These files have been abbreviated to the first 100000 lines for
    space reasons when presenting these examples.

## Standard commands to yield these output files

vdW interaction free energies:
```
julia -p 1 AARMBEM_vdW.jl mollistfilename=2graphene3p9nmby3p4nmVACconfig_periodic.txt outfilename=2graphene3p9nmby3p4nmVACvdW_periodic.out translistfilename=2graphene3p9nmby3p4nmVACtrans_periodic.txt periodicfilename=1graphene3p9nmby3p4nmlatticevecs.txt freqlistfilename=freq31.txt klistfilename=graphene3p9nmby3p4nmunitcellklist1.txt Genv=VAC
```

Thermal radiation powers:
```
julia -p 1 AARMBEM_heat.jl mollistfilename=2graphene3p9nmby3p4nmVACconfig_periodic.txt outfilename=2graphene3p9nmby3p4nmVACheat_periodic.out translistfilename=2graphene3p9nmby3p4nmVACtrans_periodic.txt periodicfilename=1graphene3p9nmby3p4nmlatticevecs.txt freqlistfilename=freq31.txt klistfilename=graphene3p9nmby3p4nmunitcellklist1.txt Genv=VAC
```

## Other possible examples for commands (not exhaustive)

Thermal radiation powers in the absence of EM retardation:
```
julia -p 1 AARMBEM_heat.jl mollistfilename=2graphene3p9nmby3p4nmVACconfig_periodic.txt outfilename=2graphene3p9nmby3p4nmVACheat_periodic.out translistfilename=2graphene3p9nmby3p4nmVACtrans_periodic.txt periodicfilename=1graphene3p9nmby3p4nmlatticevecs.txt freqlistfilename=freq31.txt klistfilename=graphene3p9nmby3p4nmunitcellklist1.txt Genv=VAC nonretarded
```
The output file name will actually be
`nonretarded_2graphene3p9nmby3p4nmVACheat_periodic.out`.

vdW interaction free energies in the absence of phonons:
```
julia -p 1 AARMBEM_vdW.jl mollistfilename=2graphene3p9nmby3p4nmVACconfig_periodic.txt outfilename=2graphene3p9nmby3p4nmVACvdW_periodic.out translistfilename=2graphene3p9nmby3p4nmVACtrans_periodic.txt periodicfilename=1graphene3p9nmby3p4nmlatticevecs.txt freqlistfilename=freq31.txt klistfilename=graphene3p9nmby3p4nmunitcellklist1.txt Genv=VAC nophonons
```
The output file name will actually be
`nophonons_2graphene3p9nmby3p4nmVACvdW_periodic.out`.

vdW interaction free energies in the absence of phonons or EM
retardation:
```
julia -p 1 AARMBEM_vdW.jl mollistfilename=2graphene3p9nmby3p4nmVACconfig_periodic.txt outfilename=2graphene3p9nmby3p4nmVACvdW_periodic.out translistfilename=2graphene3p9nmby3p4nmVACtrans_periodic.txt periodicfilename=1graphene3p9nmby3p4nmlatticevecs.txt freqlistfilename=freq31.txt klistfilename=graphene3p9nmby3p4nmunitcellklist1.txt Genv=VAC nophonons nonretarded
```
The output file name will actually be
`nonretarded_nophonons_2graphene3p9nmby3p4nmVACvdW_periodic.out`.

