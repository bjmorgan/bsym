---
title: 'bsym: A basic symmetry module"
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

# Summary

`bsym`[@bsym_version_1] is a Python module for performing basic symmetry operations on arbitrary configuration vector spaces. It consists of core classes that describe configuration vector spaces, their symmetry operations, and specific configurations of objects within these spaces. The module also contains functions for performing symmetry analyses on configurations of objects. Key functionality includes enumerating symmetry-inequivalent configurations arranged over a configuration vector space.

The primary intended use-case for `bsym` is for generating molecular or crystalline structures with fractional site occupancy, for use in modelling by chemists and solid state physicists. The code contains an efficient implementation of the algorithm described by Grau-Crespo _et al._ [@Grau-CrespoEtAl_JPhysCondensMatter2007] for enumerating the symmetry-inequivalent configurations in a disordered atomic structure. The code also includes an interface with [pymatgen](http://pymatgen.org) [@pymatgen] for performing direct symmetry analysis of crystal structures represented by pymatgen `Structure` objects.

# Acknowledgements

BJM acknowledges support from the EU FP7 program SIRBATT (contract No. 608502), from the Royal Society (UF130329), and from EPSRC grant EP/N004302/ 1.

# References

