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

`bsym` is a basic Python symmetry module. It consists of core classes that describe configuration vector spaces, their symmetry operations, and specific configurations of objects withing these spaces. The module also contains an interface for working with [`pymatgen`](http://pymatgen.org) `Structure` objects, to allow simple generation of disordered symmetry-inequivalent structures from a symmetric parent crystal structure.

The primary intended use-case for `bsym` is for generating molecular or crystalline structures with fractional site occupancy, for use in modelling by chemists and solid state physicists. The code contains an efficient implementation of the algorithm described by Grau-Crespo _et al._ [@Grau-CrespoEtAl_JPhysCondensMatter2007] for enumerating the symmetry-inequivalent configurations in a disordered atomic or crystalline structure, and an interface with [pymatgen](http://pymatgen.org) [@pymatgen] for direct symmetry analysis of pymatgen `Structure` objects.

# Acknowledgements

BJM acknowledges support from the EU FP7 program SIRBATT (contract No. 608502), from the Royal Society (UF130329), and from EPSRC grant EP/N004302/ 1.

# References

