import unittest
import numpy as np
from bsym import ColourOperation, Configuration
from unittest.mock import patch

class ColourOperationTestCase( unittest.TestCase ):
    """Tests for colour operation methods"""

    def test_symmetry_operation_is_initialised_from_a_matrix( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation( matrix, colour_mapping=mapping )
        np.testing.assert_array_equal( co.matrix, matrix )
        self.assertEqual( co.colour_mapping, mapping )

    def test_from_vector( self ):
        vector = [ 2, 3, 1 ]
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 }, { 1: 1, 0: 0 } ]
        co = ColourOperation.from_vector( vector, mapping )
        np.testing.assert_array_equal( co.matrix, np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( co.colour_mapping, mapping )

    def test_from_vector_with_label( self ):
        vector = [ 2, 3, 1 ]
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        label = 'A'
        co = ColourOperation.from_vector( vector, mapping, label=label )
        np.testing.assert_array_equal( co.matrix, np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( co.label, label )
        self.assertEqual( co.colour_mapping, mapping )

    def test_symmetry_operation_is_initialised_with_label( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        label = 'E'
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation( matrix, mapping, label=label )
        self.assertEqual( co.label, label )
        self.assertEqual( co.colour_mapping, mapping )

    def test_from_vector_counting_from_zero( self ):
        vector = [ 1, 2, 0 ]
        mapping = [ { 1: 0, 0: 1 }, { 1: 1, 0: 0 } ]
        co = ColourOperation.from_vector( vector, mapping, count_from_zero=True )
        np.testing.assert_array_equal( co.matrix, np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( co.colour_mapping, mapping )

    def test_operate_on( self ):
        matrix = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        colour_mapping = [ { 1:1, 2:2, 3:3 },
                           { 1:2, 2:3, 3:1 },
                           { 1:3, 2:2, 3:1 } ]
        co = ColourOperation( matrix, colour_mapping )
        configuration = Configuration( [ 1, 2, 3 ] )
        co.operate_on( configuration )
        np.testing.assert_array_equal( co.operate_on( configuration ).vector, np.array( [ 2, 1, 3 ] ) )

    def test_mul( self ):
        matrix_a = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        colour_mapping_a = [ { 0:1, 1:0 }, { 0:1, 1:0 } ]
        matrix_b = np.array( [ [ 0, 1 ], [ 1, 0 ] ] )
        colour_mapping_b = [ { 0:1, 1:0 }, { 0:1, 1:0 } ]
        co_a = ColourOperation( matrix_a, colour_mapping_a )
        co_b = ColourOperation( matrix_b, colour_mapping_b )
        co_c = co_a * co_b
        np.testing.assert_array_equal( co_c.matrix , np.array( [ [ 0, 1 ], [ 1, 0 ] ] ) )
        self.assertEqual( co_c.colour_mapping, [ { 0:0, 1:1 }, { 0:0, 1:1 } ] )

if __name__ == '__main__':
    unittest.main()
