import unittest
from unittest.mock import Mock, patch
from bsym import ConfigurationSpace, SymmetryGroup, SymmetryOperation, Configuration
from bsym.configuration_space import permutation_as_config_number, colourings_generator
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
    def test_unique_configurations_by_composition_calls_unique_configurations(
        self, mock_unique_permutations, mock_generate_partitions):
        """
        Test that the method calls self.unique_configurations with correct
        site_distribution for each composition.
        """
        mock_generate_partitions.return_value = [(2, 1)]
        mock_unique_permutations.return_value = iter([(2, 1), (1, 2)])
        
        config_space = ConfigurationSpace(objects=[1, 2, 3])
        with patch.object(config_space, 'unique_configurations', return_value=[]) as mock_unique_configs:
            config_space.unique_configurations_by_composition(n_species=2)
        
        self.assertEqual(mock_unique_configs.call_count, 2)
        
        # Check first call: composition (2, 1) → site_distribution {0: 2, 1: 1}
        first_call_kwargs = mock_unique_configs.call_args_list[0][1]
        self.assertEqual(first_call_kwargs['site_distribution'], {0: 2, 1: 1})
        self.assertEqual(first_call_kwargs['verbose'], False)
        self.assertEqual(first_call_kwargs['show_progress'], False)
        
        # Check second call: composition (1, 2) → site_distribution {0: 1, 1: 2}
        second_call_kwargs = mock_unique_configs.call_args_list[1][1]
        self.assertEqual(second_call_kwargs['site_distribution'], {0: 1, 1: 2})
    
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
        Test that configurations returned by unique_configurations are stored
        in the results dict under the correct composition tuple key.
        """
        mock_generate_partitions.return_value = [(2, 1)]
        mock_unique_permutations.return_value = iter([(2, 1), (1, 2)])
        
        # Create distinguishable mock configurations
        config_a = Mock(spec=Configuration)
        config_a.label = 'a'
        config_b = Mock(spec=Configuration)
        config_b.label = 'b'
        
        config_space = ConfigurationSpace(objects=[1, 2, 3])
        
        def mock_unique_configs_side_effect(site_distribution, **kwargs):
            if site_distribution == {0: 2, 1: 1}:
                return [config_a]
            elif site_distribution == {0: 1, 1: 2}:
                return [config_b]
        
        with patch.object(config_space, 'unique_configurations', 
                         side_effect=mock_unique_configs_side_effect):
            result = config_space.unique_configurations_by_composition(n_species=2)
        
        self.assertEqual(len(result), 2)
        self.assertEqual(result[(2, 1)], [config_a])
        self.assertEqual(result[(1, 2)], [config_b])
    
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
        Test that verbose=True prints per-composition and summary information.
        """
        mock_generate_partitions.return_value = [(2, 1)]
        mock_unique_permutations.return_value = iter([(2, 1), (1, 2)])
        
        mock_config = Mock(spec=Configuration)
        config_space = ConfigurationSpace(objects=[1, 2, 3])
        
        with patch.object(config_space, 'unique_configurations', return_value=[mock_config]):
            config_space.unique_configurations_by_composition(n_species=2, verbose=True)
        
        output = mock_stdout.getvalue()
        
        # Check for per-composition output
        self.assertIn('Processing composition (2, 1)', output)
        self.assertIn('Found 1 unique configurations', output)
        self.assertIn('Processing composition (1, 2)', output)
        
        # Check for summary output
        self.assertIn('Summary:', output)
        self.assertIn('Evaluated 2 compositions', output)
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
        
        # Verify tqdm was called with correct total
        mock_tqdm.assert_called_once()
        call_args = mock_tqdm.call_args
        # Should be called with total=2 (number of valid compositions)
        self.assertEqual(call_args[1]['total'], 2)  # keyword arg 'total'
        
        # Verify progress bar was updated
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
            

if __name__ == '__main__':
    unittest.main()
