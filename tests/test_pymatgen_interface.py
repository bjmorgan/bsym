import unittest
from unittest.mock import Mock, patch, call
import numpy as np
from pymatgen import Lattice, Structure
from bsym.pymatgen_interface import unique_symmetry_operations_as_vectors_from_structure, spacegroup_from_structure
from itertools import permutations
from bsym.symmetry_operation import SymmetryOperation
from bsym.spacegroup import SpaceGroup

class TestPymatgenInterface( unittest.TestCase ):

    def setUp( self ):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.0, 0.5, 0.5 ],
                             [ 0.5, 0.0, 0.5 ] ] )
        atom_list = [ 'Li' ] * len( coords )
        lattice = Lattice.from_parameters( a = 3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        self.structure = Structure( lattice, atom_list, coords )

    def test_unique_symmetry_operations_as_vectors_from_structure( self ):
        # integration test
        mappings = unique_symmetry_operations_as_vectors_from_structure( self.structure, verbose=False )
        self.assertEqual( len( mappings ), 24 )
        for l in permutations( [ 1, 2, 3, 4 ], 4 ):
            self.assertEqual( list(l) in mappings, True )

    @patch( 'bsym.pymatgen_interface.unique_symmetry_operations_as_vectors_from_structure' )
    @patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' )
    @patch( 'bsym.pymatgen_interface.SpaceGroup' )
    def test_spacegroup_from_structure( self, mock_SpaceGroup, mock_symmetry_operation_from_vector, mock_symmetry_operations_from_structure ):
        mock_symmetry_operations_from_structure.return_value=[ [ 1, 2 ], [ 2, 1 ] ]
        mock_symmetry_operation_from_vector.side_effect = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation) ]
        mock_SpaceGroup.return_value = Mock( spec=SpaceGroup )
        spacegroup = spacegroup_from_structure( self.structure )
        self.assertEqual( spacegroup, mock_SpaceGroup.return_value )
        self.assertEqual( mock_symmetry_operation_from_vector.call_args_list, [call([1, 2]), call([2, 1])] )

if __name__ == '__main__':
    unittest.main()
