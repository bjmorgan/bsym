import unittest
from unittest.mock import Mock, patch, call
from bsym import CoordinateConfigSpace, SymmetryGroup, ConfigurationSpace
import numpy as np

def mock_configuration_space_init( self, objects, symmetry_group ):
    self.objects = [ 'A', 'B' ]

class CoordinateConfigSpaceTestCase( unittest.TestCase ):

    @patch( "bsym.coordinate_config_space.ConfigurationSpace.__init__" )
    def test_coordinate_config_space_is_initialised_and_calls_super( self, mock_super_init ):
        #mock_super_init.create_autospec( ConfigurationSpace )
        coordinates = np.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ] )
        coordinate_config_space = CoordinateConfigSpace( coordinates )
        np.testing.assert_array_equal( coordinate_config_space.coordinates, coordinates )
        np.testing.assert_array_equal( mock_super_init.call_args[0][0], np.array( [ 1, 2 ] ) )
        self.assertEqual( mock_super_init.call_args[0][1], None )

    @patch( "bsym.coordinate_config_space.ConfigurationSpace.__init__", new=mock_configuration_space_init )
    def test_coordinate_config_space_is_initialised_and_sets_objects( self ):
        coordinates = np.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ] )
        coordinate_config_space = CoordinateConfigSpace( coordinates )
        np.testing.assert_array_equal( coordinate_config_space.coordinates, coordinates )
        np.testing.assert_array_equal( coordinate_config_space.objects, [ 'A', 'B' ] )

    @patch( "bsym.coordinate_config_space.ConfigurationSpace.__init__" )
    def test_coordinate_config_space_is_initialised_with_objects( self, mock_super_init ):
        objects = np.array( [ 3, 4 ] )
        coordinates = np.array( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ] )
        coordinate_config_space = CoordinateConfigSpace( coordinates, objects=objects )
        np.testing.assert_array_equal( coordinate_config_space.coordinates, coordinates )
        np.testing.assert_array_equal( mock_super_init.call_args[0][0], np.array( [ 3, 4 ] ) )
        self.assertEqual( mock_super_init.call_args[0][1], None )
    
if __name__ == '__main__':
    unittest.main()        
