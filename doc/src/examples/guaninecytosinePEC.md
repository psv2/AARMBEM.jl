# Guanine and cytosine above a PEC plane

This page links to examples of input files needed to compute vdW
interaction free energies and thermal radiation powers in a system of
compact molecular bodies. The particular system under consideration is
of guanine and cytosine at the same distance ``z`` above a PEC plane
(at ``z = 0``), and where cytosine may be rotated clockwise about the
``z``-axis through its center of mass by an angle ``\varphi``.

## Molecular configuration files

The input data files are as follows.

``B_{\mathrm{e}}`` for cytosine: [cytosineBe.txt](cytosineBe.txt)

``B_{\mathrm{I}}`` for cytosine: [cytosineBI.txt](cytosineBI.txt)

``K_{\mathrm{e}}`` for cytosine: [cytosineKe.txt](cytosineKe.txt)

``K_{\mathrm{I}}`` for cytosine: [cytosineKI_mod.txt](cytosineKI_mod.txt)

``M_{\mathrm{e}}`` for cytosine: [cytosineMe.txt](cytosineMe.txt)

``M_{\mathrm{I}}`` for cytosine: [cytosineMI.txt](cytosineMI.txt)

``Q_{\mathrm{e}}`` for cytosine: [cytosineQe.txt](cytosineQe.txt)

Atomic coordinates ``\{\vec{r}_{p}\}`` for cytosine: [cytosinexyz.txt](cytosinexyz.txt)


``B_{\mathrm{e}}`` for guanine: [guanineBe.txt](guanineBe.txt)

``B_{\mathrm{I}}`` for guanine: [guanineBI.txt](guanineBI.txt)

``K_{\mathrm{e}}`` for guanine: [guanineKe.txt](guanineKe.txt)

``K_{\mathrm{I}}`` for guanine: [guanineKI_mod.txt](guanineKI_mod.txt)

``M_{\mathrm{e}}`` for guanine: [guanineMe.txt](guanineMe.txt)

``M_{\mathrm{I}}`` for guanine: [guanineMI.txt](guanineMI.txt)

``Q_{\mathrm{e}}`` for guanine: [guanineQe.txt](guanineQe.txt)

Atomic coordinates ``\{\vec{r}_{p}\}`` for guanine: [guaninexyz.txt](guaninexyz.txt)


Overall configuration file: [guaninecytosineconfig.txt](guaninecytosineconfig.txt)

## Transformation file

Transformation file: [guaninecytosinetrans.txt](guaninecytosinetrans.txt)

## Frequency file

Frequency file: [freq31.txt](freq31.txt)

(Many more frequencies must be used to replicate the output for vdW
interaction free energies and thermal radiation powers: these can be
extracted from the corresponding output files themselves. For clarity,
in the example output files used, for vdW interactions, the input
frequencies ``\omega = \mathrm{i}|w|`` used were
``w_{n} = 10^{10 + (n - 1)/10}~\mathrm{rad/s}`` for
``n \in \{1, 2, \ldots, 81\}``, while for thermal radiation, the input
frequencies ``\omega = |w|`` used
were ``w_{n} = n \times 10^{12}~\mathrm{rad/s}`` for
``n \in \{1, 2, \ldots, 500\}``.)

## Output files

vdW interaction free energies: [guaninecytosinePECvdW.out](guaninecytosinePECvdW.out)

Thermal radiation powers: [guaninecytosinePECheat.out](guaninecytosinePECheat.out)

!!! note

    The ordering of output lines may vary from one run to another.

## Standard commands to yield these output files

vdW interaction free energies:
```
julia -p 1 AARMBEM_vdW.jl mollistfilename=guaninecytosineconfig.txt outfilename=guaninecytosinePECvdW.out translistfilename=guaninecytosinetrans.txt freqlistfilename=freq31.txt Genv=PEC
```

Thermal radiation powers:
```
julia -p 1 AARMBEM_heat.jl mollistfilename=guaninecytosineconfig.txt outfilename=guaninecytosinePECheat.out translistfilename=guaninecytosinetrans.txt freqlistfilename=freq31.txt Genv=PEC
```

## Other possible examples for commands (not exhaustive)

Thermal radiation powers in the absence of EM retardation:
```
julia -p 1 AARMBEM_heat.jl mollistfilename=guaninecytosineconfig.txt outfilename=guaninecytosinePECheat.out translistfilename=guaninecytosinetrans.txt freqlistfilename=freq31.txt Genv=PEC nonretarded
```
The output file name will actually be
`nonretarded_guaninecytosinePECheat.out`.

vdW interaction free energies in the absence of phonons:
```
julia -p 1 AARMBEM_vdW.jl mollistfilename=guaninecytosineconfig.txt outfilename=guaninecytosinePECvdW.out translistfilename=guaninecytosinetrans.txt freqlistfilename=freq31.txt Genv=PEC nophonons
```
The output file name will actually be
`nophonons_guaninecytosinePECvdW.out`.

vdW interaction free energies in the absence of phonons or EM
retardation:
```
julia -p 1 AARMBEM_vdW.jl mollistfilename=guaninecytosineconfig.txt outfilename=guaninecytosinePECvdW.out translistfilename=guaninecytosinetrans.txt freqlistfilename=freq31.txt Genv=PEC nophonons nonretarded
```
The output file name will actually be
`nonretarded_nophonons_guaninecytosinePECvdW.out`.

