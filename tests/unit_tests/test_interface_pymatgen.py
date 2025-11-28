import unittest
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.operations import SymmOp
from bsym.interface.pymatgen import (
    unique_symmetry_operations_as_vectors_from_structure, 
    space_group_from_structure, 
    parse_site_distribution, 
    unique_structure_substitutions, 
    new_structure_from_substitution, 
    configuration_space_from_structure, 
    space_group_symbol_from_structure, 
    configuration_space_from_molecule, 
    structure_cartesian_coordinates_mapping,
    molecule_cartesian_coordinates_mapping,
    unique_structure_substitutions_by_composition,
    random_unique_structure_substitutions
)

from itertools import permutations
from bsym import SymmetryOperation, Configuration, SpaceGroup, PointGroup, ConfigurationSpace

class TestPymatgenInterface(unittest.TestCase):

    def setUp(self):
        # construct a pymatgen Structure instance using the site fractional coordinates
        # face-centered cubic lattice
        coords = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.5, 0.0],
                           [0.0, 0.5, 0.5],
                           [0.5, 0.0, 0.5]])
        atom_list = ['Li'] * len(coords)
        lattice = Lattice.from_parameters(a=3.0, b=3.0, c=3.0, alpha=90, beta=90, gamma=90)
        self.structure = Structure(lattice, atom_list, coords)
        # construct a pymatgen Molecule instance
        # square molecule (D4h)
        m_coords = np.array([[0.0, 0.0, 0.0],
                              [1.0, 0.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [1.0, 1.0, 0.0]])
        molecule = Molecule(atom_list, m_coords)
        molecule = Molecule(molecule.species, molecule.cart_coords - molecule.center_of_mass)
        self.molecule = molecule 

    def test_new_structure_from_substitution(self):
        substitution_index = [2,3]
        new_species_list = ['Mg', 'Fe'] 
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

    def test_parse_site_distribution(self):
        site_distribution = {'Mg': 1, 'Li': 3}
        n, d = parse_site_distribution(site_distribution)
        for k, v in n.items():
            self.assertEqual(site_distribution[d[k]], v)

    def test_structure_cartesian_coordinates_mapping(self):
        mock_symmop = Mock(spec=SymmOp)
        new_coords = np.array([[0.5, 0.5, 0.5]])
        mock_symmop.operate_multi = Mock(return_value=new_coords)
        self.structure.lattice.get_cartesian_coords = Mock(return_value=np.array([[2.0, 2.0, 2.0]]))
        mapped_coords = structure_cartesian_coordinates_mapping(self.structure, mock_symmop)
        np.testing.assert_array_equal(mapped_coords, np.array([[2.0, 2.0, 2.0]]))
        np.testing.assert_array_equal(mock_symmop.operate_multi.call_args[0][0], self.structure.frac_coords)

    def test_molecule_cartesian_coordinates_mapping(self):
        mock_symmop = Mock(spec=SymmOp)
        new_coords = np.array([[0.5, 0.5, 0,5]])
        mock_symmop.operate_multi = Mock(return_value=new_coords)
        mapped_coords = molecule_cartesian_coordinates_mapping(self.molecule, mock_symmop) 
        np.testing.assert_array_equal(mapped_coords, new_coords)
        np.testing.assert_array_equal(mock_symmop.operate_multi.call_args[0][0], self.molecule.cart_coords)
        
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_unique_structure_substitutions_by_composition_calls_config_space_correctly(
        self, mock_new_structure, mock_config_space_from_structure):
        """Test that configuration space is created with correct parameters"""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.unique_configurations_by_composition = Mock(return_value={})
        mock_config_space_from_structure.return_value = mock_config_space
        
        unique_structure_substitutions_by_composition(
            mock_structure,
            'X',
            ['Li', 'Na'],
            atol=1e-6
        )
        
        # Verify configuration_space_from_structure called with correct args
        mock_config_space_from_structure.assert_called_once_with(
            mock_structure,
            subset=[0, 1, 2, 3],
            atol=1e-6
        )
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    def test_unique_structure_substitutions_by_composition_converts_bounds_to_indices(
        self, mock_config_space_from_structure):
        """Test that species name bounds are converted to numeric indices"""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.unique_configurations_by_composition = Mock(return_value={})
        mock_config_space_from_structure.return_value = mock_config_space
        
        bounds = {'Li': (1, 2), 'Na': (0, 2)}
        
        unique_structure_substitutions_by_composition(
            mock_structure,
            'X',
            ['Li', 'Na', 'Mg'],
            bounds=bounds
        )
        
        # Verify bounds were converted: Li->0, Na->1
        call_kwargs = mock_config_space.unique_configurations_by_composition.call_args[1]
        expected_bounds = {0: (1, 2), 1: (0, 2)}
        self.assertEqual(call_kwargs['bounds'], expected_bounds)
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    def test_unique_structure_substitutions_by_composition_raises_error_for_invalid_species_in_bounds(
        self, mock_config_space_from_structure):
        """Test ValueError raised when bounds contain species not in species_list"""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space_from_structure.return_value = mock_config_space
        
        bounds = {'K': (1, 2)}  # K not in species_list
        
        with self.assertRaises(ValueError) as context:
            unique_structure_substitutions_by_composition(
                mock_structure,
                'X',
                ['Li', 'Na'],
                bounds=bounds
            )
        
        self.assertIn("'K'", str(context.exception))
        self.assertIn("not found in species_list", str(context.exception))
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_unique_structure_substitutions_by_composition_maps_configs_to_structures(
        self, mock_new_structure, mock_config_space_from_structure):
        """Test that configurations are correctly mapped to structures with species"""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2])
        
        # Create mock configurations
        mock_config1 = Mock(spec=Configuration)
        mock_config1.tolist = Mock(return_value=[0, 1, 0])  # Li, Na, Li
        mock_config1.count = 3
        
        mock_config2 = Mock(spec=Configuration)
        mock_config2.tolist = Mock(return_value=[1, 1, 0])  # Na, Na, Li
        mock_config2.count = 2
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.unique_configurations_by_composition = Mock(
            return_value={(2, 1): [mock_config1, mock_config2]}
        )
        mock_config_space_from_structure.return_value = mock_config_space
        
        mock_structure1 = Mock(spec=Structure)
        mock_structure2 = Mock(spec=Structure)
        mock_new_structure.side_effect = [mock_structure1, mock_structure2]
        
        results = unique_structure_substitutions_by_composition(
            mock_structure,
            'X',
            ['Li', 'Na']
        )
        
        # Verify new_structure_from_substitution called with correct species
        calls = mock_new_structure.call_args_list
        self.assertEqual(len(calls), 2)
        
        # First call: [0, 1, 0] -> ['Li', 'Na', 'Li']
        self.assertEqual(calls[0][0][2], ['Li', 'Na', 'Li'])
        
        # Second call: [1, 1, 0] -> ['Na', 'Na', 'Li']
        self.assertEqual(calls[1][0][2], ['Na', 'Na', 'Li'])
        
        # Verify metadata was set
        self.assertEqual(mock_structure1.number_of_equivalent_configurations, 3)
        self.assertEqual(mock_structure2.number_of_equivalent_configurations, 2)
        
        # Verify results structure
        self.assertEqual(len(results), 1)
        self.assertIn((2, 1), results)
        self.assertEqual(results[(2, 1)], [mock_structure1, mock_structure2])
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    def test_unique_structure_substitutions_by_composition_passes_through_parameters(
        self, mock_config_space_from_structure):
        """Test that verbose and show_progress parameters are passed through"""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.unique_configurations_by_composition = Mock(return_value={})
        mock_config_space_from_structure.return_value = mock_config_space
        
        unique_structure_substitutions_by_composition(
            mock_structure,
            'X',
            ['Li', 'Na'],
            verbose=True,
            show_progress='notebook'
        )
        
        call_kwargs = mock_config_space.unique_configurations_by_composition.call_args[1]
        self.assertEqual(call_kwargs['verbose'], True)
        self.assertEqual(call_kwargs['show_progress'], 'notebook')
        
class TestRandomUniqueStructureSubstitutions(unittest.TestCase):
    """Tests for random_unique_structure_substitutions."""
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_calls_configuration_space_from_structure_correctly(
        self, mock_new_structure, mock_config_space_from_structure
    ):
        """Test that configuration space is created with correct parameters."""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[])
        mock_config_space_from_structure.return_value = mock_config_space
        
        random_unique_structure_substitutions(
            mock_structure,
            'X',
            {'Li': 2, 'Na': 2},
            n=5,
            atol=1e-6
        )
        
        mock_config_space_from_structure.assert_called_once_with(
            mock_structure,
            subset=[0, 1, 2, 3],
            atol=1e-6
        )
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_converts_site_distribution_to_numeric(
        self, mock_new_structure, mock_config_space_from_structure
    ):
        """Test that species names are converted to numeric indices."""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[])
        mock_config_space_from_structure.return_value = mock_config_space
        
        random_unique_structure_substitutions(
            mock_structure,
            'X',
            {'Li': 2, 'Na': 2},
            n=5
        )
        
        call_kwargs = mock_config_space.random_unique_configurations.call_args[1]
        self.assertEqual(call_kwargs['site_distribution'], {0: 2, 1: 2})
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_passes_parameters_to_random_unique_configurations(
        self, mock_new_structure, mock_config_space_from_structure
    ):
        """Test that n, sampling, and seed are passed through correctly."""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[])
        mock_config_space_from_structure.return_value = mock_config_space
        
        random_unique_structure_substitutions(
            mock_structure,
            'X',
            {'Li': 2, 'Na': 2},
            n=10,
            sampling='uniform',
            seed=42
        )
        
        call_kwargs = mock_config_space.random_unique_configurations.call_args[1]
        self.assertEqual(call_kwargs['n'], 10)
        self.assertEqual(call_kwargs['sampling'], 'uniform')
        self.assertEqual(call_kwargs['seed'], 42)
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_converts_configurations_to_structures(
        self, mock_new_structure, mock_config_space_from_structure
    ):
        """Test that configurations are converted to Structure objects."""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2])
        
        mock_config = Mock(spec=Configuration)
        mock_config.tolist = Mock(return_value=[0, 1, 0])  # Li, Na, Li
        mock_config.count = 3
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[mock_config])
        mock_config_space_from_structure.return_value = mock_config_space
        
        mock_new_struct = Mock(spec=Structure)
        mock_new_structure.return_value = mock_new_struct
        
        random_unique_structure_substitutions(
            mock_structure,
            'X',
            {'Li': 2, 'Na': 1},
            n=1
        )
        
        mock_new_structure.assert_called_once_with(
            mock_structure,
            [0, 1, 2],
            ['Li', 'Na', 'Li']
        )
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    @patch('bsym.interface.pymatgen.new_structure_from_substitution')
    def test_sets_number_of_equivalent_configurations(
        self, mock_new_structure, mock_config_space_from_structure
    ):
        """Test that degeneracy is set on returned structures."""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2])
        
        mock_config = Mock(spec=Configuration)
        mock_config.tolist = Mock(return_value=[0, 1, 0])
        mock_config.count = 6
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[mock_config])
        mock_config_space_from_structure.return_value = mock_config_space
        
        mock_new_struct = Mock(spec=Structure)
        mock_new_structure.return_value = mock_new_struct
        
        result = random_unique_structure_substitutions(
            mock_structure,
            'X',
            {'Li': 2, 'Na': 1},
            n=1
        )
        
        self.assertEqual(result[0].number_of_equivalent_configurations, 6)
    
    @patch('bsym.interface.pymatgen.configuration_space_from_structure')
    def test_returns_correct_number_of_structures(
        self, mock_config_space_from_structure
    ):
        """Test that the correct number of structures is returned."""
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_configs = []
        for i in range(5):
            mock_config = Mock(spec=Configuration)
            mock_config.tolist = Mock(return_value=[0, 0, 1, 1])
            mock_config.count = 2
            mock_configs.append(mock_config)
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=mock_configs)
        mock_config_space_from_structure.return_value = mock_config_space
        
        with patch('bsym.interface.pymatgen.new_structure_from_substitution') as mock_new:
            mock_new.return_value = Mock(spec=Structure)
            result = random_unique_structure_substitutions(
                mock_structure,
                'X',
                {'Li': 2, 'Na': 2},
                n=5
            )
        
        self.assertEqual(len(result), 5)
        
    def test_random_unique_structure_substitutions_exclude_file_excludes_configurations(self):
        """
        Test that configurations from exclude_file are excluded.
        """
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[])
        
        mock_excluded_config = Mock(spec=Configuration)
        mock_excluded_config.tolist = Mock(return_value=[0, 0, 1, 1])  # 4 elements to match 4 sites
        mock_excluded_configs = [mock_excluded_config]
        
        with patch('bsym.interface.pymatgen.configuration_space_from_structure', return_value=mock_config_space):
            with patch('bsym.interface.pymatgen.load_configurations', return_value=mock_excluded_configs) as mock_load:
                random_unique_structure_substitutions(
                    mock_structure,
                    'X',
                    {'Li': 2, 'Na': 2},
                    n=1,
                    exclude_file='excluded.json',
                )
                
                mock_load.assert_called_once_with('excluded.json')
                call_kwargs = mock_config_space.random_unique_configurations.call_args[1]
                self.assertEqual(call_kwargs['exclude'], mock_excluded_configs)
    
    
    def test_random_unique_structure_substitutions_exclude_file_accepts_list(self):
        """
        Test that exclude_file accepts a list of files.
        """
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[])
        
        mock_config_1 = Mock(spec=Configuration)
        mock_config_1.tolist = Mock(return_value=[0, 0, 1, 1])
        mock_configs_1 = [mock_config_1]
        
        mock_config_2 = Mock(spec=Configuration)
        mock_config_2.tolist = Mock(return_value=[0, 1, 0, 1])
        mock_configs_2 = [mock_config_2]
        
        with patch('bsym.interface.pymatgen.configuration_space_from_structure', return_value=mock_config_space):
            with patch('bsym.interface.pymatgen.load_configurations', side_effect=[mock_configs_1, mock_configs_2]) as mock_load:
                random_unique_structure_substitutions(
                    mock_structure,
                    'X',
                    {'Li': 2, 'Na': 2},
                    n=1,
                    exclude_file=['batch_1.json', 'batch_2.json'],
                )
                
                self.assertEqual(mock_load.call_count, 2)
                mock_load.assert_any_call('batch_1.json')
                mock_load.assert_any_call('batch_2.json')
                call_kwargs = mock_config_space.random_unique_configurations.call_args[1]
                self.assertEqual(call_kwargs['exclude'], mock_configs_1 + mock_configs_2)
                
    def test_random_unique_structure_substitutions_output_file_saves_configurations(self):
        """
        Test that output_file saves the generated configurations.
        """
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])
        
        mock_config = Mock(spec=Configuration)
        mock_config.tolist = Mock(return_value=[0, 0, 1, 1])
        mock_config.count = 2
        
        mock_config_space = Mock(spec=ConfigurationSpace)
        mock_config_space.random_unique_configurations = Mock(return_value=[mock_config])
        
        with patch('bsym.interface.pymatgen.configuration_space_from_structure', return_value=mock_config_space):
            with patch('bsym.interface.pymatgen.new_structure_from_substitution') as mock_new_structure:
                mock_new_structure.return_value = Mock(spec=Structure)
                with patch('bsym.interface.pymatgen.save_configurations') as mock_save:
                    random_unique_structure_substitutions(
                        mock_structure,
                        'X',
                        {'Li': 2, 'Na': 2},
                        n=1,
                        output_file='output.json',
                    )
                    
                    mock_save.assert_called_once_with([mock_config], 'output.json')
                    
    def test_random_unique_structure_substitutions_raises_on_mismatched_configuration_length(self):
        """
        Test that ValueError is raised if excluded configuration length
        doesn't match the number of sites to substitute.
        """
        mock_structure = Mock(spec=Structure)
        mock_structure.indices_from_symbol = Mock(return_value=[0, 1, 2, 3])  # 4 sites
        
        mock_config = Mock(spec=Configuration)
        mock_config.tolist = Mock(return_value=[0, 1, 0])  # 3 elements - mismatch
        
        with patch('bsym.interface.pymatgen.configuration_space_from_structure'):
            with patch('bsym.interface.pymatgen.load_configurations', return_value=[mock_config]):
                with self.assertRaises(ValueError):
                    random_unique_structure_substitutions(
                        mock_structure,
                        'X',
                        {'Li': 2, 'Na': 2},
                        n=1,
                        exclude_file='excluded.json',
                    )
        

if __name__ == '__main__':
    unittest.main()
