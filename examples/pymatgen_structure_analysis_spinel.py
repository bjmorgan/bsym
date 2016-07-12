#! /usr/bin/env python3

# Uses the interface with pymatgen to construct a 2x2x2 simple cubic lattice, then identify the symmetry inequivalent structures with 2 sites substituted.

import numpy as np
from pymatgen import Lattice, Structure
from bsym.pymatgen_interface import spacegroup_from_structure
from bsym import bsym as sym

if __name__ == '__main__':
    # construct a pymatgen Structure instance using the site fractional coordinates
    coords = np.array( [ [ 0.00, 0.00, 0.00],
                         [ 0.5, 0.5, 0.00 ],
                         [ 0.25,0.25,0.25],
                         [ 0.75,0.75,0.25],
                         [ 0.5,0.0,0.5],
                         [ 0.0,0.5,0.5],
                         [ 0.75,0.25,0.75],
                         [ 0.25,0.75,0.75] ] )
    atom_list = [ 'Li' ] * len( coords )
    lattice = Lattice.from_parameters( a = 1.0, b=1.0, c=1.0, alpha=90, beta=90, gamma=90 )         
    parent_structure = Structure( lattice, atom_list, coords )

    # generate a SpaceGroup instance with the symmetry operations for this structure
    sg = spacegroup_from_structure( parent_structure )

    # define the site occupations
    site_dist = { 1 : 1, 
                  0 : 4 }

    # find the unique configurations
    unique_configurations = sym.process.unique_configurations_from_sites( site_dist, sg, verbose=True )

    # output
    [ config.pprint() for config in unique_configurations ] 
    sitelist = sym.SiteList( coords )
    sym.process.coordinate_list_from_sitelists( configs = unique_configurations,
                                                labels = [ 0, 1 ],
                                                sitelists = [ sitelist ] )


