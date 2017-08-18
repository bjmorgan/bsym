---
title: 'bsym: A basic symmetry module'
tags:  
  - symmetry
  - configurations
  - disorder
authors:  
 - name: Benjamin J. Morgan  
   email: b.j.morgan@bath.ac.uk  
   orcid: 0000-0002-3056-8233  
   affiliation: 1  
affiliations:  
 - name: Department of Chemistry, University of Bath, Bath, BA2 7AY, United Kingdom.  
   index: 1  
date: 22 July 2017  
bibliography: paper.bib
---

A large number of problems in materials science concern the configurational disorder of atoms. Examples include describing mixing of alloys [@GanoseAndScanlon_JMaterChemA2016]; identifying prefferred arrangements in crystals of non-dilute point defects, dopants, or intercalated atoms [@MorganAndWatson_JPhysChemLett2011, @GrieshammerEtAl_PhysChemChemPhys2014, @DaltonEtAl_ChemMater2012]; or predicting crystal surface structures or the arrangements of adsorbed atoms [@MorganEtAl_JMaterChemA2016].

Computational modelling can provide useful insight into problems such as these. For a given atomic geometry, relative energies, or other such properties, can be directly calculated. Repeating these calculations across a range of possible structures can be used to identify which atom arrangements are more or less likely, or be used to construct a statistical description giving an ensemble average.
For bulk materials with approximately random disorder (ideal mixing) special quasi-random structures can be used to model the effects of disorder in a single periodic system [@ZungerEtAl_PhysRevLett1990]. Where disorder is not random, however, contributions from different atomic configurations need to be explicitly considered. One approach is to model the system energy via an effective Hamiltonian (for example, parameterised using cluster expansion methods), and then to perform Monte Carlo simulations [@LerchEtAl_ModellingSimulMaterSciEng2009, @VanDeWalleAndCeder_JPhaseEquil2002, @CASM_v0.2.1, @LudwigEtAl_PhysRevB2011]. This approach depends on the accuracy of the effective Hamiltonian, and is particularly suited to systems where short-range interactions dominate the total energy. An alternative is to explicitly consider all possible configurations of relevant atoms within a reduced configurational space (e.g. a computationally tractable supercell) [@Grau-CrespoEtAl_JPhysCondensMatter2007, @TompsettAndIslam_ChemMater2013]. 

In this approach, the overall computational cost can be greatly reduced by considering only *symmetry inequivalent* configurations. If, for each of these unique structures, the number of symmetry-equivalent configurations is known, then in addition to individual enthalpies, configurational entropies, and hence free energies, may be calculated [@Grau-CrespoEtAl_JPhysCondensMatter2007]. 

# `bsym`

`bsym` [@bsym_version_1] is a Python module for performing symmetry-based manipulations on arbitrary configuration vector spaces. 
The code includes an efficient implementation of the algorithm described by Grau-Crespo *et al.*, for enumerating symmetry inequivalent configurations [@Grau-CrespoEtAl_JPhysCondensMatter2007].
For treating the specific case of site-disorder in crystals, `bsym` includes an interface for working with `pymatgen` `Structure` objects [@OngEtAl_CompMaterSci2013], which allows simple construction of sets of symmetry operations for crystal structures, and conversion to and from a range of standard file formats for recording crystal structures and for inputs for a range of atomistic modelling codes. 
The core classes describe objects corresponding to an abstracted matrix representation of configurational vector spaces and their symmetry properties. This means `bsym` can be readily used for analysis in other classes of problems, including molecular point groups, crystal surface symmetries (2D space groups), and graph theoretical colouring problems.

# Acknowledgements

BJM acknowledges support from the EU FP7 program SIRBATT (contract No. 608502), from the Royal Society (UF130329), and from EPSRC grant EP/N004302/ 1.

# References

