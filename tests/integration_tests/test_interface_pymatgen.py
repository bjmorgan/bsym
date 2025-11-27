import unittest
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from bsym.interface.pymatgen import (
    unique_symmetry_operations_as_vectors_from_structure, 
    space_group_from_structure, 
    parse_site_distribution, 
    unique_structure_substitutions, 
    new_structure_from_substitution, 
    configuration_space_from_structure, 
    space_group_symbol_from_structure, 
    configuration_space_from_molecule,
    unique_structure_substitutions_by_composition,
    random_unique_structure_substitutions
)

from itertools import permutations
from bsym import ( SymmetryOperation, 
                   Configuration, 
                   SpaceGroup, 
                   PointGroup, 
                   ConfigurationSpace )

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
        # integration test
        mock_structure = Mock( spec=Structure )
        mock_structure.indices_from_symbol = Mock( return_value = [ 0, 1, 2 ] )
        with self.assertRaises( ValueError ):
            unique_structure_substitutions( mock_structure, 'Li', { 'A':1, 'B':1 } )
         
    def test_space_group_symbol_from_structure( self ):
        # integration test
        self.assertEqual( space_group_symbol_from_structure( self.structure ), 'Fm-3m' )
        
    def test_unique_structure_substitutions_by_composition_binary_on_square(self):
        """Test binary substitution on 4-site square lattice gives all expected compositions"""
        # Create a 4-site square lattice structure
        coords = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.0, 0.0],
                           [0.0, 0.5, 0.0],
                           [0.5, 0.5, 0.0]])
        atom_list = ['X'] * 4  # Placeholder atoms to be substituted
        lattice = Lattice.from_parameters(a=2.0, b=2.0, c=2.0, alpha=90, beta=90, gamma=90)
        parent_structure = Structure(lattice, atom_list, coords)
        
        # Perform composition-based substitution
        results = unique_structure_substitutions_by_composition(
            parent_structure,
            'X',
            ['Li', 'Na']
        )
        
        # Should have 5 compositions: (4,0), (3,1), (2,2), (1,3), (0,4)
        self.assertEqual(len(results), 5)
        
        # Check (4, 0): all Li
        self.assertIn((4, 0), results)
        self.assertEqual(len(results[(4, 0)]), 1)
        self.assertEqual(results[(4, 0)][0].composition.get_atomic_fraction('Li'), 1.0)
        self.assertEqual(results[(4, 0)][0].number_of_equivalent_configurations, 1)
        
        # Check (3, 1): 3 Li, 1 Na
        self.assertIn((3, 1), results)
        self.assertEqual(len(results[(3, 1)]), 1)
        self.assertEqual(results[(3, 1)][0].composition.get_atomic_fraction('Li'), 0.75)
        self.assertEqual(results[(3, 1)][0].number_of_equivalent_configurations, 4)
        
        # Check (2, 2): 2 Li, 2 Na - should have 2 unique structures
        self.assertIn((2, 2), results)
        self.assertEqual(len(results[(2, 2)]), 2)
        total_degeneracy_2_2 = sum(s.number_of_equivalent_configurations for s in results[(2, 2)])
        self.assertEqual(total_degeneracy_2_2, 6)  # Should be C(4,2) = 6
        # Check stoichiometry
        for s in results[(2, 2)]:
            self.assertEqual(s.composition.get_atomic_fraction('Li'), 0.5)
            self.assertEqual(s.composition.get_atomic_fraction('Na'), 0.5)
        
        # Verify one structure has adjacent arrangement, one has diagonal
        distances_squared = []
        for s in results[(2, 2)]:
            li_indices = s.indices_from_symbol('Li')
            na_indices = s.indices_from_symbol('Na')
            # Distance between the two Li atoms
            dist_sq = s.get_distance(li_indices[0], li_indices[1])**2
            distances_squared.append(dist_sq)
        distances_squared = sorted(distances_squared)
        np.testing.assert_array_almost_equal(distances_squared, [1.0, 2.0])  # adjacent=1, diagonal=√2
        
        # Check (1, 3): 1 Li, 3 Na
        self.assertIn((1, 3), results)
        self.assertEqual(len(results[(1, 3)]), 1)
        self.assertEqual(results[(1, 3)][0].composition.get_atomic_fraction('Li'), 0.25)
        self.assertEqual(results[(1, 3)][0].number_of_equivalent_configurations, 4)
        
        # Check (0, 4): all Na
        self.assertIn((0, 4), results)
        self.assertEqual(len(results[(0, 4)]), 1)
        self.assertEqual(results[(0, 4)][0].composition.get_atomic_fraction('Na'), 1.0)
        self.assertEqual(results[(0, 4)][0].number_of_equivalent_configurations, 1)
        
        # Verify total unique configurations
        total_unique = sum(len(configs) for configs in results.values())
        self.assertEqual(total_unique, 6)
        
class TestRandomUniqueStructureSubstitutions(unittest.TestCase):
    """Integration tests for random_unique_structure_substitutions."""
    
    def setUp(self):
        """Set up a simple 4x4 square lattice structure."""
        coords = np.array([[0.0, 0.0, 0.0]])
        atom_list = ['Li']
        lattice = Lattice.from_parameters(a=1.0, b=1.0, c=1.0, alpha=90, beta=90, gamma=90)
        self.parent_structure = Structure(lattice, atom_list, coords) * [4, 4, 1]
    
    def test_returns_n_unique_structures(self):
        """Test that the correct number of structures is returned."""
        result = random_unique_structure_substitutions(
            self.parent_structure,
            'Li',
            {'Na': 2, 'Li': 14},
            n=3,
            seed=42
        )
        
        self.assertEqual(len(result), 3)
    
    def test_returned_structures_are_mutually_inequivalent(self):
        """Test that returned structures are symmetry-inequivalent."""
        result = random_unique_structure_substitutions(
            self.parent_structure,
            'Li',
            {'Na': 2, 'Li': 14},
            n=5,
            seed=42
        )
        
        # Each structure should have Na at different relative positions
        na_positions = []
        for struct in result:
            na_indices = struct.indices_from_symbol('Na')
            na_coords = tuple(sorted([tuple(struct[i].frac_coords) for i in na_indices]))
            na_positions.append(na_coords)
        
        # All Na position sets should be unique
        self.assertEqual(len(na_positions), len(set(na_positions)))
    
    def test_same_seed_produces_same_results(self):
        """Test that using the same seed produces identical results."""
        result_1 = random_unique_structure_substitutions(
            self.parent_structure,
            'Li',
            {'Na': 2, 'Li': 14},
            n=3,
            seed=42
        )
        
        result_2 = random_unique_structure_substitutions(
            self.parent_structure,
            'Li',
            {'Na': 2, 'Li': 14},
            n=3,
            seed=42
        )
        
        self.assertEqual(len(result_1), len(result_2))
        for s1, s2 in zip(result_1, result_2):
            na_idx_1 = list(s1.indices_from_symbol('Na'))
            na_idx_2 = list(s2.indices_from_symbol('Na'))
            self.assertEqual(na_idx_1, na_idx_2)
    
    def test_structures_have_correct_composition(self):
        """Test that returned structures have the requested composition."""
        result = random_unique_structure_substitutions(
            self.parent_structure,
            'Li',
            {'Na': 2, 'Mg': 1, 'Li': 13},
            n=3,
            seed=42
        )
        
        for struct in result:
            composition = struct.composition.as_dict()
            self.assertEqual(composition['Na'], 2)
            self.assertEqual(composition['Mg'], 1)
            self.assertEqual(composition['Li'], 13)
    
    def test_structures_have_degeneracy_attribute(self):
        """Test that returned structures have number_of_equivalent_configurations set."""
        result = random_unique_structure_substitutions(
            self.parent_structure,
            'Li',
            {'Na': 2, 'Li': 14},
            n=3,
            seed=42
        )
        
        for struct in result:
            self.assertTrue(hasattr(struct, 'number_of_equivalent_configurations'))
            self.assertIsInstance(struct.number_of_equivalent_configurations, int)
            self.assertGreater(struct.number_of_equivalent_configurations, 0)

if __name__ == '__main__':
    unittest.main()
