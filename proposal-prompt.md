/opsx-continue (i.e. create proposal) In this project we intend to implement some methods which are described (without great detail) in scientific papers. They simulate inelastic neutron scattering (INS) measurements using molecular dynamics simulation data; rather than working directly from the trajectory data these methods use derivedphonon density-of-states (and similar) data as an intermediate representation of the system dynamics.

The first method is described in https://pubs.acs.org/jctcce/article/16/12/7702/617412/Simulation-of-Inelastic-Neutron-Scattering-Spectra and uses the phonon DOS in place of atomic displacement tensors <u> in a fully-isotropic approximation.

The second method is described in https://www.nature.com/articles/s41598-021-86771-5 and goes further to include cross-correlation terms in the energy-dependent input.

The first capability will be to generate the atom-projected density-of-states from trajectory data; external libraries could be used as long as they are well-maintained and fit well in a modern Python ecosystem. This can be validated against the pDOS from harmonic calculations computed with Euphonic from a consistent force calculator such as the ASE Lennard-Jones implementation.

From there we intend to implement the INS intensity calculations described in the above-mentioned papers and benchmark against existing codes.