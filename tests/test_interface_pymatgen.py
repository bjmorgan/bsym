import unittest
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np
from pymatgen import Lattice, Structure, Molecule
from bsym.interface.pymatgen import unique_symmetry_operations_as_vectors_from_structure, space_group_from_structure, parse_site_distribution, unique_structure_substitutions, new_structure_from_substitution, configuration_space_from_structure, space_group_symbol_from_structure, configuration_space_from_molecule

from itertools import permutations
from bsym import SymmetryOperation, Configuration, SpaceGroup, PointGroup, ConfigurationSpace

class TestPymatgenInterface( unittest.TestCase ):

    def setUp( self ):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                             [ 0.5, 0.5, 0.0 ],
                             [ 0.0, 0.5, 0.5 ],
                             [ 0.5, 0.0, 0.5 ] ] )
        atom_list = [ 'Li' ] * len( coords )
        lattice = Lattice.from_parameters( a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90 )
        self.structure = Structure( lattice, atom_list, coords )
        # construct a pymatgen Molecule instance
        # square molecule (D4h)
        m_coords = np.array( [ [ 0.0, 0.0, 0.0 ],
                               [ 1.0, 0.0, 0.0 ],
                               [ 0.0, 1.0, 0.0 ],
                               [ 1.0, 1.0, 0.0 ] ] )
        molecule = Molecule( atom_list, m_coords )
        molecule = Molecule( molecule.species, molecule.cart_coords - molecule.center_of_mass )
        self.molecule = molecule 

    def test_unique_symmetry_operations_as_vectors_from_structure( self ):
        # integration test
        mappings = unique_symmetry_operations_as_vectors_from_structure( self.structure, verbose=False )
        self.assertEqual( len( mappings ), 24 )
        for l in permutations( [ 1, 2, 3, 4 ], 4 ):
            self.assertEqual( list(l) in mappings, True )

    def test_unique_symmetry_operations_as_vectors_from_structure_hex( self ):
        # integration test
        coords = np.array( [ [ 0.666667, 0.333334, 0.498928 ],
                             [ 0.333334, 0.666667, 0.998928 ],
                             [ 0.666667, 0.333334, 0.876081 ],
                             [ 0.333334, 0.666667, 0.376081 ] ] )
        atom_list = [ 'Zn' ] * 2 + [ 'O' ] * 2
        lattice = Lattice.from_parameters( a=2.0, b=2.0, c=3.265986324, alpha=90, beta=90, gamma=120 )
        structure = Structure( lattice, atom_list, coords )
        mappings = unique_symmetry_operations_as_vectors_from_structure( structure, verbose=False )
        self.assertEqual( len( mappings ), 2 )
        for l in [ [ 1, 2, 3, 4 ], [ 2, 1, 4, 3 ] ]:
            self.assertEqual( l in mappings, True )

    @patch( 'bsym.interface.pymatgen.unique_symmetry_operations_as_vectors_from_structure' )
    @patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' )
    @patch( 'bsym.interface.pymatgen.SpaceGroup' )
    def test_space_group_from_structure( self, mock_SpaceGroup, mock_symmetry_operation_from_vector, mock_symmetry_operations_from_structure ):
        mock_symmetry_operations_from_structure.return_value=[ [ 1, 2 ], [ 2, 1 ] ]
        mock_symmetry_operation_from_vector.side_effect = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation) ]
        mock_SpaceGroup.return_value = Mock( spec=SpaceGroup )
        space_group = space_group_from_structure( self.structure )
        self.assertEqual( space_group, mock_SpaceGroup.return_value )
        self.assertEqual( mock_symmetry_operation_from_vector.call_args_list, [call([1, 2]), call([2, 1])] )

    @patch( 'bsym.interface.pymatgen.space_group_from_structure' )
    @patch( 'bsym.interface.pymatgen.ConfigurationSpace' )
    def test_configuration_space_from_structure( self, mock_ConfigurationSpace, mock_space_group_from_structure ):
        mock_space_group = Mock( spec=SpaceGroup )
        mock_space_group_from_structure.return_value = mock_space_group
        mock_configspace = Mock( spec=ConfigurationSpace )
        mock_ConfigurationSpace.return_value = mock_configspace
        config_space = configuration_space_from_structure( self.structure )
        self.assertEqual( config_space, mock_configspace )
        mock_space_group_from_structure.assert_called_with( self.structure, subset=None, atol=1e-5 )
        mock_ConfigurationSpace.assert_called_with( objects=[ 1, 2, 3, 4 ], symmetry_group=mock_space_group )

    @patch( 'bsym.interface.pymatgen.point_group_from_molecule' )
    @patch( 'bsym.interface.pymatgen.ConfigurationSpace' )
    def test_configuration_space_from_molecule( self, mock_ConfigurationSpace, mock_point_group_from_molecule ):
        mock_point_group = Mock( spec=PointGroup )
        mock_point_group_from_molecule.return_value = mock_point_group
        mock_configspace = Mock( spec=ConfigurationSpace )
        mock_ConfigurationSpace.return_value = mock_configspace
        config_space = configuration_space_from_molecule( self.molecule )
        self.assertEqual( config_space, mock_configspace )
        mock_point_group_from_molecule.assert_called_with( self.molecule, subset=None, atol=1e-5 )
        mock_ConfigurationSpace.assert_called_with( objects=[ 1, 2, 3, 4 ], symmetry_group=mock_point_group )
 
    @patch( 'bsym.interface.pymatgen.unique_symmetry_operations_as_vectors_from_structure' )
    @patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' )
    @patch( 'bsym.interface.pymatgen.SpaceGroup' )
    def test_space_group_from_structure_with_subset_calls_with_subset( self, mock_SpaceGroup, mock_symmetry_operation_from_vector, mock_symmetry_operations_from_structure ):
        mock_symmetry_operations_from_structure.return_value=[ [ 1, 2, ], [ 2, 1 ] ]
        mock_symmetry_operation_from_vector.side_effect = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation) ]
        mock_SpaceGroup.return_value = Mock( spec=SpaceGroup )
        subset = [ 0 ]
        atol = 1e-5
        space_group = space_group_from_structure( self.structure, subset=subset )
        mock_symmetry_operations_from_structure.assert_called_once_with( self.structure, subset=subset, atol=atol )
 
    def test_unique_structure_colourings( self ):
        # integration test
        c = configuration_space_from_molecule( self.molecule )
        uc = c.unique_colourings( [ 0, 1 ] )
        self.assertEqual( len( uc ), 6 )
 
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
        np.testing.assert_array_equal( np.array( sorted( [ s.full_configuration_degeneracy for s in ns ] ) ), np.array( [ 1, 2, 4, 4, 4 ] ) )

    def test_unique_structure_substitutions_in_two_steps_gives_full_degeneracies( self ):
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
        us = unique_structure_substitutions( parent_structure, 'Li', { 'Na':1, 'Li':15 } )
        ns = unique_structure_substitutions( us[0], 'Li', { 'Mg':1, 'Li':14 } )
        self.assertEqual( len( ns ), 5 )
        distances = np.array( sorted( [ s.get_distance( s.indices_from_symbol('Mg')[0], s.indices_from_symbol('Na')[0] ) for s in ns ] ) )
        np.testing.assert_array_almost_equal( distances, np.array( [ 0.75    ,  1.06066 ,  1.5     ,  1.677051,  2.12132 ] ) )
        np.testing.assert_array_equal( np.array( sorted( [ s.number_of_equivalent_configurations for s in ns ] ) ), np.array( [ 1, 2, 4, 4, 4 ] ) )
        np.testing.assert_array_equal( np.array( sorted( [ s.full_configuration_degeneracy for s in ns ] ) ), np.array( [ 16, 32, 64, 64, 64 ] ) )

    def test_unique_structure_substitutions_with_mismatched_site_distribution_raises_ValueError( self ):
        mock_structure = Mock( spec=Structure )
        mock_structure.indices_from_symbol = Mock( return_value = [ 0, 1, 2 ] )
        with self.assertRaises( ValueError ):
            unique_structure_substitutions( mock_structure, 'Li', { 'A':1, 'B':1 } )
         
    def test_pymatgen_structure_can_be_patched( self ):
        with self.assertRaises( AttributeError ):
            self.structure.number_of_equivalent_configurations
        with self.assertRaises( AttributeError ):
            self.structure.full_configuration_degeneracy
 
    def test_new_structure_from_substitution( self ):
        substitution_index = [ 2,3 ]
        new_species_list = [ 'Mg', 'Fe' ] 
        s_new = new_structure_from_substitution( self.structure, substitution_index, new_species_list ) 
        self.assertEqual( s_new[2].species_string, 'Mg' )
        self.assertEqual( s_new[3].species_string, 'Fe' )

    def test_new_structure_from_substitution_raises_ValueError_with_oversize_index( self ):
        substitution_index = [ 0, 1, 2, 3, 4 ]
        new_species_list = [ 'Mg', 'Fe' ]
        with self.assertRaises( ValueError ):
            new_structure_from_substitution( self.structure, substitution_index, new_species_list )

    def test_new_structure_from_substitution_raises_ValueError_with_invalid_index( self ):
        substitution_index = [ 2, 4 ]
        new_species_list = [ 'Mg', 'Fe' ]
        with self.assertRaises( ValueError ):
            new_structure_from_substitution( self.structure, substitution_index, new_species_list )

    def test_space_group_symbol_from_structure( self ):
        self.assertEqual( space_group_symbol_from_structure( self.structure ), 'Fm-3m' )

class TestPymatgenAPI( unittest.TestCase ):
    
    def test_parse_site_distribution( self ):
        sd = { 'Li': 2, 'Mg': 4 }
        sd_numeric, sd_mapping = parse_site_distribution( sd )    
        for k, v in sd_numeric.items():
            self.assertEqual( sd[ sd_mapping[ k ] ], v )

if __name__ == '__main__':
    unittest.main()
