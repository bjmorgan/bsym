import unittest
from unittest.mock import Mock, patch
from bsym import ConfigurationSpace, SymmetryGroup, SymmetryOperation, Configuration
from bsym.configuration_space import permutation_as_config_number, colourings_generator
import numpy as np

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

class ConfigurationSpaceModuleFunctionsTestCase( unittest.TestCase ):
      
    def test_permutation_as_config_number( self ):
        self.assertEqual( permutation_as_config_number( [ 1, 1, 0, 0, 1 ] ), 11001 )

    def test_colourings_generator( self ):
        colourings = list( colourings_generator( [ 1, 0 ], dim=3 ) )
        expected_colourings = [ [1, 1, 1], 
                                [0, 1, 1], [1, 0, 1], [1, 1, 0], 
                                [0, 0, 1], [0, 1, 0], [1, 0, 0], 
                                [0, 0, 0] ]
        for c in colourings:
            self.assertEqual( c in expected_colourings, True )
        for ec in expected_colourings:
            self.assertEqual( ec in colourings, True )
        
if __name__ == '__main__':
    unittest.main()
