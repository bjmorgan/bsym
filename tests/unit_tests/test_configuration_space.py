import unittest
from unittest.mock import Mock, patch
from bsym import ConfigurationSpace,SymmetryOperation, Configuration
from bsym.configuration_space import permutation_as_config_number, colourings_generator, apply_species_mapping
from bsym.symmetry_group import SymmetryGroup
import numpy as np
import io

class ConfigurationSpaceTestCase( unittest.TestCase ):

    def test_configuration_space_is_initialised( self ):
        mock_symmetry_group = Mock( spec=SymmetryGroup )
        mock_symmetry_operations = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation ) ]
        mock_symmetry_operations[0].matrix = np.array( np.zeros( (3,3) ) )
        mock_symmetry_operations[1].matrix = np.array( np.zeros( (3,3) ) )
        mock_symmetry_group.symmetry_operations = mock_symmetry_operations
        object_list = [ 'A', 'B', 'C' ]
        configuration_space = ConfigurationSpace( symmetry_group=mock_symmetry_group, objects=object_list )
        self.assertEqual( configuration_space.symmetry_group, mock_symmetry_group )
        self.assertEqual( configuration_space.objects, object_list )

    def test_configuration_space_initialisation_raises_valueerror_if_dimensions_are_inconsistent( self ):
        mock_symmetry_group = Mock( spec=SymmetryGroup )
        mock_symmetry_operations = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation ) ]
        mock_symmetry_operations[0].matrix = np.array( np.zeros( (3,3) ) )
        mock_symmetry_operations[1].matrix = np.array( np.zeros( (3,3) ) )
        mock_symmetry_group.symmetry_operations = mock_symmetry_operations
        object_list = [ 'A', 'B' ]
        with self.assertRaises( ValueError ):
            ConfigurationSpace( symmetry_group=mock_symmetry_group, objects=object_list ) 

    def test_configuration_space_initialised_with_no_symmetry_group( self ):
        object_list = [ 'A', 'B' ]
        configuration_space = ConfigurationSpace( objects=object_list )
        self.assertEqual( configuration_space.symmetry_group.size, 1 )
        self.assertEqual( configuration_space.symmetry_group.symmetry_operations[0].label, 'E' )
        np.testing.assert_array_equal( configuration_space.symmetry_group.symmetry_operations[0].matrix, np.array( [[1,0],[0,1]] ) )

    def test_configuration_space_initialised_with_no_symmetry_group_creates_sym_op_with_ints( self ):
        object_list = [ 'A', 'B' ]
        configuration_space = ConfigurationSpace( objects=object_list )
        self.assertEqual( issubclass( configuration_space.symmetry_group.symmetry_operations[0].matrix.dtype.type, np.integer ), True )

    def test_unique_configurations( self ):
        object_list = [ 1, 1, 2 ]
        configuration_space = ConfigurationSpace( objects=object_list )
        site_distribution = { 1:2, 2:1 }
        mock_configuration = Mock( spec=Configuration )
        configuration_space.enumerate_configurations = Mock( return_value=[ mock_configuration ] )
        with patch( 'bsym.configuration_space.flatten_list' ) as mock_flatten_list:
            mock_flatten_list.return_value = [ 1, 1, 2 ] 
            with patch( 'bsym.configuration_space.unique_permutations' ) as mock_unique_permutations:
                mock_unique_permutations.return_value = [ [ 1, 1, 2 ], [ 1, 2, 1 ], [ 2, 1, 1 ] ]
                configurations = configuration_space.unique_configurations( site_distribution )
                mock_unique_permutations.assert_called_with( [ 1, 1, 2 ] )
            mock_flatten_list.assert_called_with( [ [1, 1 ], [ 2 ] ] )
        configuration_space.enumerate_configurations.assert_called_with(
                          mock_unique_permutations(), verbose=False )
        self.assertEqual( configurations, [ mock_configuration ] )

    def test_unique_colourings( self ):
        object_list = [ 1, 1, 2 ]
        configuration_space = ConfigurationSpace( objects=object_list )
        mock_configuration = Mock( spec=Configuration )
        mock_configuration.dim = 3
        configuration_space.enumerate_configurations = Mock( return_value=[ mock_configuration ] )
        with patch( 'bsym.configuration_space.colourings_generator' ) as mock_colourings_generator:
            mock_colourings_generator.return_values = [ [ 1, 1, 2 ], [ 1, 2, 1 ], [ 2, 1, 1 ] ]
            colourings = configuration_space.unique_colourings( colours=[ 1, 2 ] )
            mock_colourings_generator.assert_called_with( [1, 2], mock_configuration.dim )    
        configuration_space.enumerate_configurations.assert_called_with(
                          mock_colourings_generator(), verbose=False )
        self.assertEqual( colourings, [ mock_configuration ] )
            
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_returns_dict(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test that method returns a dict with composition tuples as keys
        and lists of Configuration objects as values.
        """
        mock_generate_partitions.return_value = [(2, 0)]
        mock_unique_permutations.return_value = iter([(2, 0)])
        mock_config = Mock(spec=Configuration)
        
        config_space = ConfigurationSpace(objects=[1, 2])
        with patch.object(config_space, 'unique_configurations', return_value=[mock_config]):
            result = config_space.unique_configurations_by_composition(n_species=2)
        
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 1)
        self.assertIn((2, 0), result)
        self.assertIsInstance(result[(2, 0)], list)
        self.assertEqual(result[(2, 0)], [mock_config])
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_excludes_zero_counts_from_site_distribution(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test that species with count=0 are excluded from site_distribution
        passed to unique_configurations.
        """
        mock_generate_partitions.return_value = [(3, 0, 1)]
        mock_unique_permutations.return_value = iter([(3, 0, 1)])
        
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            config_space.unique_configurations_by_composition(n_species=3)
        
        call_kwargs = mock_unique_configs.call_args[1]
        # Species 1 (index 1) has count 0, so should be excluded
        self.assertEqual(call_kwargs['site_distribution'], {0: 3, 2: 1})
        self.assertNotIn(1, call_kwargs['site_distribution'])
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_stores_configs_under_correct_keys(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test that configurations are stored under correct composition keys.
        Canonical composition is analyzed, non-canonical are relabeled.
        """
        mock_generate_partitions.return_value = [(2, 1)]
        mock_unique_permutations.return_value = iter([(2, 1), (1, 2)])
        
        # Create mock configuration for canonical
        config_a = Mock(spec=Configuration)
        config_a.label = 'a'
        config_a.vector = np.array([0, 0, 1])  # 2 of species 0, 1 of species 1
        config_a.count = 3
        
        config_space = ConfigurationSpace(objects=[1, 2, 3])
        
        def mock_unique_configs_side_effect(site_distribution, **kwargs):
            if site_distribution == {0: 2, 1: 1}:  # canonical (2, 1)
                return [config_a]
            else:
                raise ValueError(f"Unexpected site_distribution: {site_distribution}")
        
        with patch.object(config_space, 'unique_configurations', 
                          side_effect=mock_unique_configs_side_effect):
            result = config_space.unique_configurations_by_composition(n_species=2)
        
        # Both compositions should be in results
        self.assertEqual(len(result), 2)
        self.assertIn((2, 1), result)
        self.assertIn((1, 2), result)
        
        # Canonical (2, 1) should have the original mock
        self.assertEqual(result[(2, 1)], [config_a])
        
        # Non-canonical (1, 2) should have relabeled version
        self.assertEqual(len(result[(1, 2)]), 1)
        relabeled = result[(1, 2)][0]
        
        # Check relabelled config has correct species mapping: [0,0,1] → [1,1,0]
        np.testing.assert_array_equal(relabeled.vector, np.array([1, 1, 0]))
        self.assertEqual(relabeled.count, 3)  # Count should be preserved
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    @patch('bsym.configuration_space.satisfies_bounds')
    def test_unique_configurations_by_composition_filters_by_bounds(
        self, mock_satisfies_bounds, mock_unique_permutations, mock_generate_partitions):
        """
        Test that compositions are filtered by bounds before processing.
        """
        mock_generate_partitions.return_value = [(3, 1), (2, 2)]
        # First partition generates 2 permutations, second generates 1
        mock_unique_permutations.side_effect = [
            iter([(3, 1), (1, 3)]),  # permutations of (3, 1)
            iter([(2, 2)])           # permutations of (2, 2)
        ]
        
        # Only allow (3, 1) and (2, 2) to pass bounds check
        def satisfies_bounds_side_effect(composition_dict, bounds):
            if composition_dict == {0: 3, 1: 1}:
                return True
            elif composition_dict == {0: 1, 1: 3}:
                return False
            elif composition_dict == {0: 2, 1: 2}:
                return True
            return False
        
        mock_satisfies_bounds.side_effect = satisfies_bounds_side_effect
        
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        bounds = {0: (2, 3), 1: (1, 2)}
        
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            result = config_space.unique_configurations_by_composition(n_species=2, bounds=bounds)
        
        # Should only process 2 compositions (3,1) and (2,2)
        self.assertEqual(mock_unique_configs.call_count, 2)
        self.assertEqual(len(result), 2)
        self.assertIn((3, 1), result)
        self.assertIn((2, 2), result)
        self.assertNotIn((1, 3), result)
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    @patch('bsym.configuration_space.satisfies_bounds')
    def test_unique_configurations_by_composition_no_bounds_skips_filtering(
        self, mock_satisfies_bounds, mock_unique_permutations, mock_generate_partitions):
        """
        Test that when bounds=None, satisfies_bounds is not called.
        """
        mock_generate_partitions.return_value = [(2, 2)]
        mock_unique_permutations.return_value = iter([(2, 2)])
        
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        with patch.object(config_space, 'unique_configurations', return_value=[]):
            config_space.unique_configurations_by_composition(n_species=2, bounds=None)
        
        mock_satisfies_bounds.assert_not_called()
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_single_species(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test edge case: n_species=1 (only one possible composition).
        """
        mock_generate_partitions.return_value = [(4,)]
        mock_unique_permutations.return_value = iter([(4,)])
        
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            result = config_space.unique_configurations_by_composition(n_species=1)
        
        self.assertEqual(len(result), 1)
        self.assertIn((4,), result)
        
        # Check that site_distribution is correct
        call_kwargs = mock_unique_configs.call_args[1]
        self.assertEqual(call_kwargs['site_distribution'], {0: 4})
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_more_species_than_sites(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test edge case: n_species > n_sites (some species must have count 0).
        """
        mock_generate_partitions.return_value = [(2, 0, 0), (1, 1, 0)]
        mock_unique_permutations.side_effect = [
            iter([(2, 0, 0), (0, 2, 0), (0, 0, 2)]),
            iter([(1, 1, 0), (1, 0, 1), (0, 1, 1)])
        ]
        
        config_space = ConfigurationSpace(objects=[1, 2])
        
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            result = config_space.unique_configurations_by_composition(n_species=3)
        
        self.assertEqual(len(result), 6)
        
        # Check that zero counts are excluded from site_distributions
        for call in mock_unique_configs.call_args_list:
            site_dist = call[1]['site_distribution']
            # No species should have count 0 in site_distribution
            for count in site_dist.values():
                self.assertNotEqual(count, 0)
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_unique_configurations_by_composition_verbose_false(
        self, mock_stdout, mock_unique_permutations, mock_generate_partitions):
        """
        Test that verbose=False doesn't print output.
        """
        mock_generate_partitions.return_value = [(2, 0)]
        mock_unique_permutations.return_value = iter([(2, 0)])
        
        config_space = ConfigurationSpace(objects=[1, 2])
        
        with patch.object(config_space, 'unique_configurations', return_value=[]):
            config_space.unique_configurations_by_composition(n_species=2, verbose=False)
        
        output = mock_stdout.getvalue()
        self.assertEqual(output, '')
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_unique_configurations_by_composition_verbose_true(
        self, mock_stdout, mock_unique_permutations, mock_generate_partitions):
        """
        Test that verbose=True prints per-partition and summary information.
        """
        mock_generate_partitions.return_value = [(2, 1)]
        mock_unique_permutations.return_value = iter([(2, 1), (1, 2)])
        
        mock_config = Mock(spec=Configuration)
        mock_config.vector = np.array([0, 0, 1])
        mock_config.count = 1
        
        config_space = ConfigurationSpace(objects=[1, 2, 3])
        
        with patch.object(config_space, 'unique_configurations', return_value=[mock_config]):
            config_space.unique_configurations_by_composition(n_species=2, verbose=True)
        
        output = mock_stdout.getvalue()
        
        # Check for per-partition output
        self.assertIn('Processing partition (2, 1)', output)
        self.assertIn('Found 1 unique configurations', output)
        
        # Check for summary output
        self.assertIn('Summary:', output)
        self.assertIn('Analyzed 1 partitions', output)
        self.assertIn('Generated 2 compositions', output)
        self.assertIn('Total unique configurations: 2', output)
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    @patch('bsym.configuration_space.tqdm')
    def test_unique_configurations_by_composition_show_progress_true(
        self, mock_tqdm, mock_unique_permutations, mock_generate_partitions):
        """
        Test that show_progress=True creates progress bar.
        """
        mock_generate_partitions.return_value = [(2, 1)]
        mock_unique_permutations.return_value = iter([(2, 1), (1, 2)])
        
        # Create mock progress bar
        mock_progress_bar = Mock()
        mock_tqdm.return_value = mock_progress_bar
        
        config_space = ConfigurationSpace(objects=[1, 2, 3])
        
        with patch.object(config_space, 'unique_configurations', return_value=[]):
            config_space.unique_configurations_by_composition(n_species=2, show_progress=True)
        
        # Verify tqdm was called (without checking total since it's now indeterminate)
        mock_tqdm.assert_called_once()
        
        # Verify progress bar was updated (once per composition: (2,1) and (1,2))
        self.assertEqual(mock_progress_bar.update.call_count, 2)
        
        # Verify progress bar was closed
        mock_progress_bar.close.assert_called_once()
    
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_passes_show_progress_to_unique_configurations(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test that show_progress parameter is passed through to 
        self.unique_configurations for nested progress bars.
        """
        mock_generate_partitions.return_value = [(2, 0)]
        mock_unique_permutations.side_effect = lambda x: iter([(2, 0)])
        
        config_space = ConfigurationSpace(objects=[1, 2])
        
        # Test with show_progress=True
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            config_space.unique_configurations_by_composition(n_species=2, show_progress=True)
        
        call_kwargs = mock_unique_configs.call_args[1]
        self.assertEqual(call_kwargs['show_progress'], True)
        
        # Test with show_progress=False
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            config_space.unique_configurations_by_composition(n_species=2, show_progress=False)
        
        call_kwargs = mock_unique_configs.call_args[1]
        self.assertEqual(call_kwargs['show_progress'], False)
        
    @patch('bsym.configuration_space.generate_partitions')
    @patch('bsym.configuration_space.unique_permutations')
    def test_unique_configurations_by_composition_only_analyzes_canonical_compositions(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test that symmetry analysis is only performed on canonical compositions.
        
        For a 2-site system with 2 species:
        - Compositions: (2, 0), (1, 1), (0, 2)
        - (2, 0) and (0, 2) are related by species exchange
        - Only (2, 0) and (1, 1) should undergo symmetry analysis
        - All three compositions should be present in the output
        """
        mock_generate_partitions.return_value = [(2, 0), (1, 1)]
        mock_unique_permutations.side_effect = [
            iter([(2, 0), (0, 2)]),  # permutations of partition (2, 0)
            iter([(1, 1)])           # permutations of partition (1, 1)
        ]
        
        # Create mock configurations
        mock_config_2_0 = Mock(spec=Configuration)
        mock_config_2_0.vector = np.array([0, 0])  # Both sites are species 0
        mock_config_2_0.count = 1
        
        mock_config_1_1 = Mock(spec=Configuration)
        mock_config_1_1.vector = np.array([0, 1])  # One of each species
        mock_config_1_1.count = 1
        
        mock_config_0_2 = Mock(spec=Configuration)  # Shouldn't be needed if working correctly
        mock_config_0_2.vector = np.array([1, 1])  # Both sites are species 1
        mock_config_0_2.count = 1
        
        config_space = ConfigurationSpace(objects=[1, 2])
        
        # Mock unique_configurations to return different configs for different compositions
        def mock_unique_configs_side_effect(site_distribution, **kwargs):
            if site_distribution == {0: 2}:  # composition (2, 0)
                return [mock_config_2_0]
            elif site_distribution == {0: 1, 1: 1}:  # composition (1, 1)
                return [mock_config_1_1]
            elif site_distribution == {1: 2}:  # composition (0, 2) - shouldn't be called
                return [mock_config_0_2]
            else:
                return []  # Default fallback
        
        with patch.object(config_space, 'unique_configurations',
                         side_effect=mock_unique_configs_side_effect) as mock_unique_configs:
            result = config_space.unique_configurations_by_composition(n_species=2)
        
        # Verify unique_configurations called exactly twice (not for (0, 2))
        self.assertEqual(mock_unique_configs.call_count, 2)
        
        # Verify calls for (2, 0) and (1, 1)
        call_args_list = mock_unique_configs.call_args_list
        called_site_distributions = [call[1]['site_distribution'] for call in call_args_list]
        self.assertIn({0: 2}, called_site_distributions)
        self.assertIn({0: 1, 1: 1}, called_site_distributions)
        
        # Verify (0, 2) was NOT analyzed (no site_distribution={1: 2})
        self.assertNotIn({1: 2}, called_site_distributions)
        
        # Verify all three compositions present in output
        self.assertEqual(len(result), 3)
        self.assertIn((2, 0), result)
        self.assertIn((1, 1), result)
        self.assertIn((0, 2), result)
        
        # Verify each composition returns a list
        self.assertIsInstance(result[(2, 0)], list)
        self.assertIsInstance(result[(1, 1)], list)
        self.assertIsInstance(result[(0, 2)], list)
        
        # Verify (0, 2) has configurations (generated from (2, 0))
        self.assertEqual(len(result[(0, 2)]), len(result[(2, 0)]))
    

class ConfigurationSpaceModuleFunctionsTestCase( unittest.TestCase ):
      
    def test_permutation_as_config_number(self):
        self.assertEqual( permutation_as_config_number([1, 1, 0, 0, 1]), 11001)

    def test_colourings_generator(self):
        colourings = list(colourings_generator([1, 0], dim=3))
        expected_colourings = [(1, 1, 1),
                            (0, 1, 1), (1, 0, 1), (1, 1, 0),
                            (0, 0, 1), (0, 1, 0), (1, 0, 0),
                            (0, 0, 0)]
        for c in colourings:
            self.assertEqual(c in expected_colourings, True)
        for ec in expected_colourings:
            self.assertEqual(ec in colourings, True)
            
    def test_apply_species_mapping_identity(self):
        """Identity mapping should return equivalent configuration."""
        config = Configuration([0, 1, 0, 1])
        config.count = 3
        mapping = [0, 1]  # No change
        
        result = apply_species_mapping(config, mapping)
        
        np.testing.assert_array_equal(result.vector, np.array([0, 1, 0, 1]))
        self.assertEqual(result.count, 3)
        self.assertIsNot(result, config)  # Should be a new object
    
    def test_apply_species_mapping_binary_swap(self):
        """Swap two species: 0↔1."""
        config = Configuration([0, 0, 1])
        config.count = 2
        mapping = [1, 0]  # Swap species 0 and 1
        
        result = apply_species_mapping(config, mapping)
        
        np.testing.assert_array_equal(result.vector, np.array([1, 1, 0]))
        self.assertEqual(result.count, 2)
    
    def test_apply_species_mapping_ternary_rotation(self):
        """Three-species rotation: 0→1, 1→2, 2→0."""
        config = Configuration([0, 1, 2, 0])
        config.count = 5
        mapping = [1, 2, 0]  # Rotate: 0→1, 1→2, 2→0
        
        result = apply_species_mapping(config, mapping)
        
        np.testing.assert_array_equal(result.vector, np.array([1, 2, 0, 1]))
        self.assertEqual(result.count, 5)
    
    def test_apply_species_mapping_partial_swap(self):
        """Swap first two of three species: 0↔1, keep 2."""
        config = Configuration([0, 1, 2, 1, 0])
        config.count = 10
        mapping = [1, 0, 2]  # Swap 0↔1, keep 2
        
        result = apply_species_mapping(config, mapping)
        
        np.testing.assert_array_equal(result.vector, np.array([1, 0, 2, 0, 1]))
        self.assertEqual(result.count, 10)
    
    def test_apply_species_mapping_preserves_count(self):
        """Degeneracy count should always be preserved."""
        config = Configuration([0, 1])
        config.count = 42
        mapping = [1, 0]
        
        result = apply_species_mapping(config, mapping)
        
        self.assertEqual(result.count, 42)
        
    def test_enumerate_configurations_finds_unique_configs(self):
        """Test that enumerate_configurations correctly identifies unique configurations."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])  # Swap first two
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        config_space = ConfigurationSpace(objects=[1, 2, 3], symmetry_group=sg)
        
        # [1,0,0] and [0,1,0] are equivalent under s1, should get 1 unique config
        # [0,0,1] is different
        permutations = iter([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
        
        unique_configs = config_space.enumerate_configurations(permutations)
        
        self.assertEqual(len(unique_configs), 2)  # Two unique configs
        
    def test_enumerate_configurations_sets_correct_counts(self):
        """Test that degeneracy counts are set correctly."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        config_space = ConfigurationSpace(objects=[1, 2, 3], symmetry_group=sg)
        
        permutations = iter([(1, 0, 0)])
        unique_configs = config_space.enumerate_configurations(permutations)
        
        # [1,0,0] has 2 equivalents: [1,0,0] and [0,1,0]
        self.assertEqual(unique_configs[0].count, 2)

class TestConfigurationSpaceRandomUniqueConfigurations(unittest.TestCase):
    """Tests for ConfigurationSpace.random_unique_configurations and helpers."""
    
    def test_random_unique_configurations_returns_n_configurations(self):
        """
        Test that random_unique_configurations returns n configurations
        when n unique configurations are found.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config_1 = Mock(spec=Configuration)
        mock_config_1.as_bytes.return_value = b'config1'
        mock_config_1.get_byte_equivalents.return_value = {b'config1'}
        
        mock_config_2 = Mock(spec=Configuration)
        mock_config_2.as_bytes.return_value = b'config2'
        mock_config_2.get_byte_equivalents.return_value = {b'config2'}
        
        with patch.object(config_space, '_generate_random_configuration',
                          side_effect=[mock_config_1, mock_config_2]):
            result = config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=2,
                sampling='degeneracy_weighted',
            )
        
        self.assertEqual(len(result), 2)
    
    def test_random_unique_configurations_skips_equivalent_configurations(self):
        """
        Test that configurations equivalent to already-seen configurations
        are skipped.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config_a = Mock(spec=Configuration)
        mock_config_a.as_bytes.return_value = b'config_a'
        mock_config_a.get_byte_equivalents.return_value = {b'config_a', b'config_b'}
        
        mock_config_b = Mock(spec=Configuration)
        mock_config_b.as_bytes.return_value = b'config_b'  # Equivalent to config_a
        
        mock_config_c = Mock(spec=Configuration)
        mock_config_c.as_bytes.return_value = b'config_c'
        mock_config_c.get_byte_equivalents.return_value = {b'config_c'}
        
        with patch.object(config_space, '_generate_random_configuration',
                          side_effect=[mock_config_a, mock_config_b, mock_config_c]):
            result = config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=2,
                sampling='degeneracy_weighted',
            )
        
        self.assertEqual(len(result), 2)
        self.assertIn(mock_config_a, result)
        self.assertIn(mock_config_c, result)
        self.assertNotIn(mock_config_b, result)
    
    def test_random_unique_configurations_sets_count_attribute(self):
        """
        Test that returned configurations have their count attribute
        set to the number of equivalent configurations.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config = Mock(spec=Configuration)
        mock_config.as_bytes.return_value = b'config'
        mock_config.get_byte_equivalents.return_value = {b'equiv_1', b'equiv_2'}
        
        with patch.object(config_space, '_generate_random_configuration',
                          return_value=mock_config):
            result = config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=1,
                sampling='degeneracy_weighted',
            )
        
        self.assertEqual(result[0].count, 2)
    
    def test_random_unique_configurations_raises_for_invalid_sampling(self):
        """
        Test that an invalid sampling value raises ValueError.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        with self.assertRaises(ValueError):
            config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=1,
                sampling='invalid_option',
            )
    
    def test_random_unique_configurations_passes_seed_to_random_generator(self):
        """
        Test that the seed parameter is used to initialise the random generator.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config = Mock(spec=Configuration)
        mock_config.as_bytes.return_value = b'config'
        mock_config.get_byte_equivalents.return_value = {b'config'}
        
        with patch('bsym.configuration_space.np.random.default_rng') as mock_rng_constructor:
            mock_rng = Mock()
            mock_rng_constructor.return_value = mock_rng
            
            with patch.object(config_space, '_generate_random_configuration',
                              return_value=mock_config):
                config_space.random_unique_configurations(
                    site_distribution={1: 2, 0: 2},
                    n=1,
                    sampling='degeneracy_weighted',
                    seed=42,
                )
            
            mock_rng_constructor.assert_called_once_with(42)
    
    def test_random_unique_configurations_uniform_rejects_based_on_degeneracy(self):
        """
        Test that uniform sampling rejects configurations with probability
        proportional to their degeneracy.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config = Mock(spec=Configuration)
        mock_config.as_bytes.return_value = b'config'
        mock_config.get_byte_equivalents.return_value = {b'equiv_1', b'equiv_2'}  # degeneracy = 2
        
        mock_rng = Mock()
        mock_rng.random.return_value = 0.3  # < 0.5 (1/degeneracy = 1/2), so should accept
        
        with patch('bsym.configuration_space.np.random.default_rng', return_value=mock_rng):
            with patch.object(config_space, '_generate_random_configuration',
                              return_value=mock_config):
                result = config_space.random_unique_configurations(
                    site_distribution={1: 2, 0: 2},
                    n=1,
                    sampling='uniform',
                )
        
        self.assertEqual(len(result), 1)
    
    def test_random_unique_configurations_uniform_rejection_does_not_add_to_seen(self):
        """
        Test that when uniform sampling rejects a configuration,
        it is not added to the seen set (can be found again later).
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config = Mock(spec=Configuration)
        mock_config.as_bytes.return_value = b'config'
        mock_config.get_byte_equivalents.return_value = {b'equiv_1', b'equiv_2'}  # degeneracy = 2
        
        mock_rng = Mock()
        mock_rng.random.side_effect = [0.7, 0.3]  # First reject, then accept
        
        with patch('bsym.configuration_space.np.random.default_rng', return_value=mock_rng):
            with patch.object(config_space, '_generate_random_configuration',
                            return_value=mock_config) as mock_generate:
                result = config_space.random_unique_configurations(
                    site_distribution={1: 2, 0: 2},
                    n=1,
                    sampling='uniform',
                )
            
                self.assertEqual(mock_generate.call_count, 2)
        
        self.assertEqual(len(result), 1)
    
    def test_generate_random_configuration_returns_configuration_with_correct_distribution(self):
        """
        Test that _generate_random_configuration returns a Configuration
        with the correct count of each species.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        mock_rng = Mock()
        
        with patch('bsym.configuration_space._select_random_indices',
                   return_value=np.array([0, 1])):
            result = config_space._generate_random_configuration(
                site_distribution={1: 2, 0: 2},
                rng=mock_rng,
            )
        
        self.assertIsInstance(result, Configuration)
        result_list = result.tolist()
        self.assertEqual(result_list.count(0), 2)
        self.assertEqual(result_list.count(1), 2)
    
    
    def test_generate_random_configuration_uses_select_random_indices(self):
        """
        Test that _generate_random_configuration uses _select_random_indices
        to select positions for each species.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        mock_rng = Mock()
        
        with patch('bsym.configuration_space._select_random_indices',
                   return_value=np.array([1, 3])) as mock_select:
            result = config_space._generate_random_configuration(
                site_distribution={1: 2, 0: 2},
                rng=mock_rng,
            )
        
        mock_select.assert_called_once()
        self.assertEqual(result.tolist(), [0, 1, 0, 1])
    
    
    def test_generate_random_configuration_passes_rng_to_select_random_indices(self):
        """
        Test that _generate_random_configuration passes the rng
        to _select_random_indices.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        mock_rng = Mock()
        
        with patch('bsym.configuration_space._select_random_indices',
                return_value=np.array([0, 1])) as mock_select:
            config_space._generate_random_configuration(
                site_distribution={1: 2, 0: 2},
                rng=mock_rng,
            )
        
        args, _ = mock_select.call_args
        self.assertIs(args[2], mock_rng)
        
    def test_random_unique_configurations_excludes_provided_configurations(self):
        """
        Test that configurations in the exclude list are not returned.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config_a = Mock(spec=Configuration)
        mock_config_a.as_bytes.return_value = b'config_a'
        mock_config_a.get_byte_equivalents.return_value = {b'config_a'}
        
        mock_config_b = Mock(spec=Configuration)
        mock_config_b.as_bytes.return_value = b'config_b'
        mock_config_b.get_byte_equivalents.return_value = {b'config_b'}
        
        mock_config_c = Mock(spec=Configuration)
        mock_config_c.as_bytes.return_value = b'config_c'
        mock_config_c.get_byte_equivalents.return_value = {b'config_c'}
        
        # Exclude config_a
        mock_excluded = Mock(spec=Configuration)
        mock_excluded.get_byte_equivalents.return_value = {b'config_a'}
        
        with patch.object(config_space, '_generate_random_configuration',
                        side_effect=[mock_config_a, mock_config_b, mock_config_c]):
            result = config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=2,
                sampling='degeneracy_weighted',
                exclude=[mock_excluded],
            )
        
        self.assertEqual(len(result), 2)
        self.assertIn(mock_config_b, result)
        self.assertIn(mock_config_c, result)
        self.assertNotIn(mock_config_a, result)
        
    def test_random_unique_configurations_excludes_equivalents_of_excluded_configurations(self):
        """
        Test that configurations equivalent to those in the exclude list
        are also not returned.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        # config_a and config_b are equivalent (both in exclude's equivalents)
        mock_config_a = Mock(spec=Configuration)
        mock_config_a.as_bytes.return_value = b'config_a'
        
        mock_config_b = Mock(spec=Configuration)
        mock_config_b.as_bytes.return_value = b'config_b'
        
        mock_config_c = Mock(spec=Configuration)
        mock_config_c.as_bytes.return_value = b'config_c'
        mock_config_c.get_byte_equivalents.return_value = {b'config_c'}
        
        mock_config_d = Mock(spec=Configuration)
        mock_config_d.as_bytes.return_value = b'config_d'
        mock_config_d.get_byte_equivalents.return_value = {b'config_d'}
        
        # Excluded configuration has both config_a and config_b as equivalents
        mock_excluded = Mock(spec=Configuration)
        mock_excluded.get_byte_equivalents.return_value = {b'config_a', b'config_b'}
        
        with patch.object(config_space, '_generate_random_configuration',
                        side_effect=[mock_config_a, mock_config_b, mock_config_c, mock_config_d]):
            result = config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=2,
                sampling='degeneracy_weighted',
                exclude=[mock_excluded],
            )
        
        self.assertEqual(len(result), 2)
        self.assertIn(mock_config_c, result)
        self.assertIn(mock_config_d, result)
        self.assertNotIn(mock_config_a, result)
        self.assertNotIn(mock_config_b, result)
        
    def test_random_unique_configurations_with_exclude_none_works_as_default(self):
        """
        Test that exclude=None behaves the same as not passing exclude.
        """
        config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
        
        mock_config = Mock(spec=Configuration)
        mock_config.as_bytes.return_value = b'config'
        mock_config.get_byte_equivalents.return_value = {b'config'}
        
        with patch.object(config_space, '_generate_random_configuration',
                        return_value=mock_config):
            result = config_space.random_unique_configurations(
                site_distribution={1: 2, 0: 2},
                n=1,
                sampling='degeneracy_weighted',
                exclude=None,
            )
        
        self.assertEqual(len(result), 1)
            

if __name__ == '__main__':
    unittest.main()
