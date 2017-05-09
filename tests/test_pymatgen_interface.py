import unittest
import numpy as np
from pymatgen import Lattice, Structure
from bsym.pymatgen_interface import unique_symmetry_operations_as_vectors_from_structure
from itertools import permutations

class TestPymatgenInterface( unittest.TestCase ):

    def test_unique_symmetry_operations_as_vectors_from_structure( self ):
        # integration test
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.0, 0.5, 0.5 ],
                             [ 0.5, 0.0, 0.5 ] ] )
        atom_list = [ 'Li' ] * len( coords )
        lattice = Lattice.from_parameters( a = 3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        structure = Structure( lattice, atom_list, coords )
        mappings = unique_symmetry_operations_as_vectors_from_structure( structure, verbose=False )
        self.assertEqual( len( mappings ), 24 )
        for l in permutations( [ 1, 2, 3, 4 ], 4 ):
            self.assertEqual( list(l) in mappings, True )

if __name__ == '__main__':
    unittest.main()
