import unittest
import numpy as np
from bsym import ColourOperation
from unittest.mock import patch

class ColourOperationTestCase( unittest.TestCase ):
    """Tests for colour operation methods"""

    def test_symmetry_operation_is_initialised_from_a_matrix( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation( matrix, colour_mapping=mapping )
        np.testing.assert_array_equal( co.matrix, matrix )
        self.assertEqual( co.colour_mapping, mapping )

    def test_from_vector( self ):
        vector = [ 2, 3, 1 ]
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation.from_vector( vector, mapping )
        np.testing.assert_array_equal( co.matrix, np.matrix( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( co.colour_mapping, mapping )

    def test_from_vector_with_label( self ):
        vector = [ 2, 3, 1 ]
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        label = 'A'
        co = ColourOperation.from_vector( vector, mapping, label=label )
        np.testing.assert_array_equal( co.matrix, np.matrix( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( co.label, label )
        self.assertEqual( co.colour_mapping, mapping )

    def test_symmetry_operation_is_initialised_with_label( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        label = 'E'
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation( matrix, mapping, label=label )
        self.assertEqual( co.label, label )
        self.assertEqual( co.colour_mapping, mapping )

    def test_from_vector_counting_from_zero( self ):
        vector = [ 1, 2, 0 ]
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation.from_vector( vector, mapping, count_from_zero=True )
        np.testing.assert_array_equal( co.matrix, np.matrix( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( co.colour_mapping, mapping )

if __name__ == '__main__':
    unittest.main()
