import unittest
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np
from pymatgen import Lattice, Structure
from bsym.interface.pymatgen import unique_symmetry_operations_as_vectors_from_structure, spacegroup_from_structure, parse_site_distribution, unique_structure_substitutions
from itertools import permutations
from bsym import SymmetryOperation, Configuration, SpaceGroup, ConfigurationSpace

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

    @patch( 'bsym.interface.pymatgen.unique_symmetry_operations_as_vectors_from_structure' )
    @patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' )
    @patch( 'bsym.interface.pymatgen.SpaceGroup' )
    def test_spacegroup_from_structure( self, mock_SpaceGroup, mock_symmetry_operation_from_vector, mock_symmetry_operations_from_structure ):
        mock_symmetry_operations_from_structure.return_value=[ [ 1, 2 ], [ 2, 1 ] ]
        mock_symmetry_operation_from_vector.side_effect = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation) ]
        mock_SpaceGroup.return_value = Mock( spec=SpaceGroup )
        spacegroup = spacegroup_from_structure( self.structure )
        self.assertEqual( spacegroup, mock_SpaceGroup.return_value )
        self.assertEqual( mock_symmetry_operation_from_vector.call_args_list, [call([1, 2]), call([2, 1])] )

    @patch( 'bsym.interface.pymatgen.unique_symmetry_operations_as_vectors_from_structure' )
    @patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' )
    @patch( 'bsym.interface.pymatgen.SpaceGroup' )
    def test_spacegroup_from_structure_with_subset_calls_with_subset( self, mock_SpaceGroup, mock_symmetry_operation_from_vector, mock_symmetry_operations_from_structure ):
        mock_symmetry_operations_from_structure.return_value=[ [ 1, 2, ], [ 2, 1 ] ]
        mock_symmetry_operation_from_vector.side_effect = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation) ]
        mock_SpaceGroup.return_value = Mock( spec=SpaceGroup )
        subset = [ 0 ]
        spacegroup = spacegroup_from_structure( self.structure, subset=subset )
        mock_symmetry_operations_from_structure.assert_called_once_with( self.structure, subset=subset )
  
    def test_unique_structure_substitutions( self ):
        # integration test
        # Create a pymatgen structure with 16 sites in a 4x4 square grid
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.25, 0.0, 0.0 ],
                             [ 0.5, 0., 0.0 ],
                             [ 0.75, 0.0, 0.0 ],
                             [ 0.0, 0.25, 0.0 ],
                             [ 0.25, 0.25, 0.0 ],
                             [ 0.5, 0.25, 0.0 ],
                             [ 0.75, 0.25, 0.0 ],
                             [ 0.0, 0.5, 0.0 ],
                             [ 0.25, 0.5, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.75, 0.5, 0.0 ],
                             [ 0.0, 0.75, 0.0 ],
                             [ 0.25, 0.75, 0.0 ],
                             [ 0.5, 0.75, 0.0 ],
                             [ 0.75, 0.75, 0.0 ] ] )
        atom_list = [ 'Li' ] * len( coords )
        lattice = Lattice.from_parameters( a = 3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        parent_structure = Structure( lattice, atom_list, coords )
        parent_structure.replace( 0, 'O' ) # substitute one site with 'O'
        ns = unique_structure_substitutions( parent_structure, 'Li', { 'Na':1, 'Li':14 } )
        self.assertEqual( len( ns ), 5 )
        distances = np.array( sorted( [ s.get_distance( s.indices_from_symbol('O')[0], s.indices_from_symbol('Na')[0] ) for s in ns ] ) )
        np.testing.assert_array_almost_equal( distances, np.array( [ 0.75    ,  1.06066 ,  1.5     ,  1.677051,  2.12132 ] ) )
        np.testing.assert_array_equal( np.array( sorted( [ s.number_of_equivalent_configurations for s in ns ] ) ), np.array( [ 1, 2, 4, 4, 4 ] ) )
 
    def test_parse_site_distribution( self ):
        sd = { 'Li': 2, 'Mg': 4 }
        sd_numeric, sd_mapping = parse_site_distribution( sd )    
        for k, v in sd_numeric.items():
            self.assertEqual( sd[ sd_mapping[ k ] ], v )
    
    def test_pymatgen_structure_can_be_patched( self ):
        with self.assertRaises( AttributeError ):
            self.structure.number_of_equivalent_configurations
     
if __name__ == '__main__':
    unittest.main()
