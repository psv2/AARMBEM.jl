# Historical background

The RMB framework, of which AARMBEM.jl is a computational
implementation, was originally developed as an extension of the
many-body dispersion (MBD) framework for computing van der Waals (vdW)
interaction energies. The details of this page can be found in work by
Venkataram et al [^VenkataramARXIV2020]; this page simply summarizes
the key points.

The MBD framework was developed to compute ground state vdW
interaction energies in atoms and in molecular systems of many
different sizes, combining ab-initio accounts of local chemical
changes to electron densities in bonded atoms relative to free atoms
[^TkatchenkoPRL2009] with Coulomb potentials coupling each atom to
account for long-range interactions [^TkatchenkoPRL2012]. It has been
applied to small dimers, clusters of silicon or carbon, carbon
allotropes, biological molecules, organic molecular crystals, and many
other systems in order to determine binding energies, structural
properties, vibrational spectra, and other properties.

The RMB framework initially [^VenkataramPRL2017] extended this
framework to account for EM retardation and the possibility that the
atom-scale bodies are interacting not in vacuum but in the presence of
a macroscopic body which can be modeled with a continuum
susceptibility, and yielded predictions of vdW interaction energies in
the presence of a gold substrate or a gold cone. Later developments in
the RMB framework accounted for the contributions of phonons to EM
response in atomistic structures, yielding novel predictions of
thermal radiation powers [^VenkataramPRL2018] and vdW interaction free
energies [^VenkataramSCIADV2019] of molecules and low-dimensional
media in the presence of a PEC plane; however, for reasons that this
and other pages make clear, this came at the cost of lesser generality
in treating continuous media with respect to the practical
implementation, though nothing about the theory had changed.

## Nomenclature of "atomism" versus "continuum"

The notions of "atomistic" or "molecular" bodies versus "macroscopic"
or "continuous" bodies are not absolute, especially given that
ultimately, all bodies are made of atoms. As a rule of thumb, bodies
that are smaller than about 5 nanometers in at least one dimension or
feature size, or those that are less than 1 nanometer apart from each
other at the smallest separation, must be treated in an ab-initio
manner incorporating atom-scale effects, so we call such bodies
"atomistic" or "molecular"; this includes even low-dimensional
materials like carbon nanotubes, carbyne wires, or graphene sheets of
infinite extent. If none of the above conditions hold, then local bulk
continuum susceptibility models may suffice, and we call such bodies
"macroscopic" or "continuous". We use the terms "mesoscopic" and
"mesoscale" to refer to situations where there are feature sizes and
separations on the order of 1-100 nanometers in which certain bodies
must be treated using ab-initio atom-scale methods while others may be
treated accurately using local bulk material response models.

## Extension to general macroscopic bodies

When phonons are included, the current implementation of the RMB
framework in AARMBEM.jl can only consider interactions of atom-scale
systems in vacuum or in the presence of a single macroscopic body,
namely a PEC plane. However, nothing about the RMB framework
inherently forbids consideration of more complicated macroscopic
systems, and a promising candidate for extending the AARMBEM.jl code
to account for this may come through combination with the
[SCUFF-EM](https://homerreid.github.io/scuff-em-documentation/) API,
though there may be other candidates too with respect to computational
solvers for macroscopic EM.

[^VenkataramARXIV2020]: Prashanth S. Venkataram, Jan Hermann, Alexandre Tkatchenko, and Alejandro W. Rodriguez. "Fluctuational Electrodynamics in Atomic and Macroscopic Systems: van der Waals Interactions and Radiative Heat Transfer". [arXiv:2005.04083](https://arxiv.org/abs/2005.04083)

[^TkatchenkoPRL2009]: Alexandre Tkatchenko and Matthias Scheffler. "Accurate Molecular Van Der Waals Interactions from Ground-State Electron Density and Free-Atom Reference Data". [Phys. Rev. Lett. **102**, 073005 (2009)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.102.073005).

[^TkatchenkoPRL2012]: Alexandre Tkatchenko, Robert A. DiStasio, Jr., Roberto Car, and Matthias Scheffler. "Accurate and Efficient Method for Many-Body van der Waals Interactions". [Phys. Rev. Lett. **108**, 236402 (2012)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.236402).

[^VenkataramPRL2017]: Prashanth S. Venkataram, Jan Hermann, Alexandre Tkatchenko, and Alejandro W. Rodriguez. "Unifying Microscopic and Continuum Treatments of van der Waals and Casimir Interactions". [Phys. Rev. Lett. **118**, 266802 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.266802).

[^VenkataramPRL2018]: Prashanth S. Venkataram, Jan Hermann, Alexandre Tkatchenko, and Alejandro W. Rodriguez. "Phonon-Polariton Mediated Thermal Radiation and Heat Transfer among Molecules and Macroscopic Bodies: Nonlocal Electromagnetic Response at Mesoscopic Scales". [Phys. Rev. Lett. **121**, 045901 (2018)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.045901).

[^VenkataramSCIADV2019]: Prashanth S. Venkataram, Jan Hermann, Teerit J. Vongkovit, Alexandre Tkatchenko, and Alejandro W. Rodriguez. "Impact of nuclear vibrations on van der Waals and Casimir interactions at zero and finite temperature". [*Sci. Adv.* **5**, 11, eaaw0456 (2019)](https://advances.sciencemag.org/content/5/11/eaaw0456).