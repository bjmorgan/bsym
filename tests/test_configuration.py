import unittest
from unittest.mock import Mock, patch
from bsym.configuration import Configuration
from bsym import SymmetryOperation
import numpy as np

class TestConfiguration( unittest.TestCase ):

    def setUp( self ):
        self.configuration = Configuration( [ 1, 0, 0 ] )

    def test_matches_returns_true_for_a_match( self ):
        self.assertEqual( self.configuration.matches( self.configuration ), True )

    def test_matches_returns_false_for_a_non_match( self ):
        other_configuration = Configuration( [0, 0, 1] )
        self.assertEqual( self.configuration.matches( other_configuration ), False )

    def test_matches_raises_TypeError_for_invalid_type( self ):
        other_configuration = 'foo'
        with self.assertRaises( TypeError ):
            self.configuration.matches( other_configuration )

    def test_is_equivalent_to_if_equivalent( self ):
        test_configuration = Configuration( [0, 1, 0] )
        symmetry_operations = [ Mock( spec=SymmetryOperation ) ]
        symmetry_operations[0].operate_on = Mock( return_value=test_configuration )
        self.assertEqual( self.configuration.is_equivalent_to( test_configuration, symmetry_operations ), True )

    def test_is_equivalent_to_if_not_equivalent( self ):
        test_configuration = Configuration( [ 0, 1, 0 ] )
        symmetry_operations = [ Mock( spec=SymmetryOperation ) ]
        symmetry_operations[0].operate_on = Mock( return_value=Configuration( [ 0, 0, 1 ] ) )
        self.assertEqual( self.configuration.is_equivalent_to( test_configuration, symmetry_operations ), False )

    def test_is_in_list( self ):
        configuration_list = [ Configuration( [ 0, 1, 0 ] ), 
                               Configuration( [ 1, 0, 0 ] ) ]
        self.configuration.matches = Mock( return_value=True )
        self.assertEqual( self.configuration.is_in_list( configuration_list ), True )

    def test_is_in_list_fails( self ):
        configuration_list = [ Configuration( [ 0, 1, 0 ] ), 
                               Configuration( [ 1, 0, 0 ] ) ]
        self.configuration.matches = Mock( return_value=False )
        self.assertEqual( self.configuration.is_in_list( configuration_list ), False )

    def test_has_equivalent_in_list( self ):
        configuration_list = [ Configuration( [ 0, 1, 0 ] ), 
                               Configuration( [ 1, 0, 0 ] ) ]
        symmetry_operations = [ Mock( spec=SymmetryOperation ) ]
        self.configuration.is_equivalent_to = Mock( return_value=True )
        self.assertEqual( self.configuration.has_equivalent_in_list( configuration_list, symmetry_operations ), True )

    def test_has_equivalent_in_list_fails( self ):
        configuration_list = [ Configuration( [ 0, 1, 0 ] ), 
                               Configuration( [ 1, 0, 0 ] ) ]
        symmetry_operations = [ Mock( spec=SymmetryOperation ) ]
        self.configuration.is_equivalent_to = Mock( return_value=False )
        self.assertEqual( self.configuration.has_equivalent_in_list( configuration_list, symmetry_operations ), False )

    def test_set_lowest_numeric_representation( self ):
        symmetry_operations = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation ) ]
        c1, c2 = Mock( spec=Configuration ), Mock( spec=Configuration )
        c1.as_number = 4
        c2.as_number = 2
        symmetry_operations[0].operate_on = Mock( return_value = c1 )
        symmetry_operations[1].operate_on = Mock( return_value = c2 )
        self.configuration.set_lowest_numeric_representation( symmetry_operations )
        self.assertEqual( self.configuration.lowest_numeric_representation, 2 )

    def test_numeric_equivalents( self ):
        symmetry_operations = [ Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation ) ]
        c1, c2 = Mock( spec=Configuration ), Mock( spec=Configuration )
        c1.as_number = 4
        c2.as_number = 2
        symmetry_operations[0].operate_on = Mock( return_value = c1 )
        symmetry_operations[1].operate_on = Mock( return_value = c2 )
        self.assertEqual( self.configuration.numeric_equivalents( symmetry_operations ), [ 4, 2 ] )

    def test_as_number( self ):
        with patch( 'bsym.configuration.Configuration.tolist' ) as mock_tolist:
            mock_tolist.side_effect = [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]
            self.assertEqual( Configuration( [ 1, 0, 0 ] ).as_number, 100 )
            self.assertEqual( Configuration( [ 0, 1, 0 ] ).as_number, 10 )
            self.assertEqual( Configuration( [ 0, 0, 1 ] ).as_number, 1 )

    def test_from_tuple( self ):
        np.testing.assert_array_equal( Configuration.from_tuple( ( 1, 1, 0 ) ).vector, 
                                       Configuration( [ 1, 1, 0 ] ).vector )

    def test_tolist( self ):
        self.assertEqual( self.configuration.tolist(), [ 1, 0, 0 ] )

    def test_position( self ):
        self.assertEqual( self.configuration.position( 0 ), [ 1, 2 ] )

    def test_map_objects( self ):
        self.assertEqual( self.configuration.map_objects( [ 'A', 'B', 'C' ] ), { 1: [ 'A' ], 0: [ 'B', 'C' ] } )  

    def test_map_objects_with_incompatible_object_list_raises_ValueError( self ):
        with self.assertRaises( ValueError ):
            self.configuration.map_objects( [ 'A', 'B' ] )

if __name__ == '__main__':
    unittest.main()
