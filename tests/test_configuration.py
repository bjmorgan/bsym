import unittest
from unittest.mock import Mock
from bsym.configuration import Configuration
from bsym.symmetry_operation import SymmetryOperation

class TestConfiguration( unittest.TestCase ):

    def setUp( self ):
        self.configuration = Configuration( [ 1, 0, 0 ] )

    def test_matches_returns_true_for_a_match( self ):
        self.assertEqual( self.configuration.matches( self.configuration ), True )

    def test_matches_returns_false_for_a_non_match( self ):
        other_configuration = Configuration( [ 0, 0, 1 ] )
        self.assertEqual( self.configuration.matches( other_configuration ), False )

    def test_is_equivalent_to_if_equivalent( self ):
        test_configuration = Configuration( [ 0, 1, 0 ] )
        symmetry_operations = [ Mock( spec=SymmetryOperation ) ]
        symmetry_operations[0].operate_on = Mock( return_value=test_configuration )
        self.assertEqual( self.configuration.is_equivalent_to( test_configuration, symmetry_operations ), True )

    def test_is_equivalent_to_if_not_equivalent( self ):
        test_configuration = Configuration( [ 0, 1, 0 ] )
        symmetry_operations = [ Mock( spec=SymmetryOperation ) ]
        symmetry_operations[0].operate_on = Mock( return_value=Configuration( [ 0, 0, 1 ] ) )
        self.assertEqual( self.configuration.is_equivalent_to( test_configuration, symmetry_operations ), False )

    def test_is_in_list( self ):
        configuration_list = [ Configuration( [ 0, 1, 0 ] ), Configuration( [ 1, 0, 0 ] ) ]
        self.configuration.matches = Mock( return_value=True )
        self.assertEqual( self.configuration.is_in_list( configuration_list ), True )

    def test_is_in_list_fails( self ):
        configuration_list = [ Configuration( [ 0, 1, 0 ] ), Configuration( [ 0, 0, 1 ] ) ]
        self.configuration.matches = Mock( return_value=False )
        self.assertEqual( self.configuration.is_in_list( configuration_list ), False )

if __name__ == '__main__':
    unittest.main()
