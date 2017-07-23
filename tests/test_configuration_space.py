import unittest
from unittest.mock import Mock
from bsym import ConfigurationSpace, SymmetryGroup, SymmetryOperation
from bsym.configuration_space import permutation_as_config_number
import numpy as np

class ConfigurationSpaceTestCase( unittest.TestCase ):

    def test_configuration_space_is_initialised( self ):
        mock_symmetry_group = Mock( spec=SymmetryGroup )
        mock_symmetry_operations = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation ) ]
        mock_symmetry_operations[0].matrix = np.matrix( np.zeros( (3,3) ) )
        mock_symmetry_operations[1].matrix = np.matrix( np.zeros( (3,3) ) )
        mock_symmetry_group.symmetry_operations = mock_symmetry_operations
        object_list = [ 'A', 'B', 'C' ]
        configuration_space = ConfigurationSpace( symmetry_group=mock_symmetry_group, objects=object_list )
        self.assertEqual( configuration_space.symmetry_group, mock_symmetry_group )
        self.assertEqual( configuration_space.objects, object_list )

    def test_configuration_space_initialisation_raises_valueerror_if_dimensions_are_inconsistent( self ):
        mock_symmetry_group = Mock( spec=SymmetryGroup )
        mock_symmetry_operations = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation ) ]
        mock_symmetry_operations[0].matrix = np.matrix( np.zeros( (3,3) ) )
        mock_symmetry_operations[1].matrix = np.matrix( np.zeros( (3,3) ) )
        mock_symmetry_group.symmetry_operations = mock_symmetry_operations
        object_list = [ 'A', 'B' ]
        with self.assertRaises( ValueError ):
            ConfigurationSpace( symmetry_group=mock_symmetry_group, objects=object_list ) 

    def test_configuration_space_initialised_with_no_symmetry_group( self ):
        object_list = [ 'A', 'B' ]
        configuration_space = ConfigurationSpace( objects=object_list )
        self.assertEqual( configuration_space.symmetry_group.size, 1 )
        self.assertEqual( configuration_space.symmetry_group.symmetry_operations[0].label, 'E' )
        np.testing.assert_array_equal( configuration_space.symmetry_group.symmetry_operations[0].matrix, np.matrix( [[1,0],[0,1]] ) )

    def test_configuration_space_initialised_with_no_symmetry_group_creates_sym_op_with_ints( self ):
        object_list = [ 'A', 'B' ]
        configuration_space = ConfigurationSpace( objects=object_list )
        self.assertEqual( issubclass( configuration_space.symmetry_group.symmetry_operations[0].matrix.dtype.type, np.integer ), True )

class ConfigurationSpaceModuleFunctionsTestCase( unittest.TestCase ):
      
    def test_permutation_as_config_number( self ):
        self.assertEqual( permutation_as_config_number( [ 1, 1, 0, 0, 1 ] ), 11001 )

if __name__ == '__main__':
    unittest.main()
