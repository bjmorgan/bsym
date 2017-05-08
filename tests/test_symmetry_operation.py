import unittest
import numpy as np
from bsym.symmetry_operation import SymmetryOperation

class SymmetryOperationTestCase( unittest.TestCase ):
    """Tests for symmetry operation functions"""

    def test_symmetry_operation_is_initialised_from_a_matrix( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        np.testing.assert_array_equal( so.matrix, matrix )

    def test_symmetry_operation_is_initialised_from_an_array( self ):
        array = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( array )
        np.testing.assert_array_equal( so.matrix, np.matrix( array ) )

    def test_symmetry_operation_is_initialised_from_a_list( self ):
        this_list = [ [ 1, 0 ], [ 0, 1 ] ]
        so = SymmetryOperation( this_list )
        np.testing.assert_array_equal( so.matrix, np.matrix( this_list ) )

    def test_symmetry_operation_raises_typeerror_for_invalid_type( self ):
        objects = [ 'foo', 1, None ]
        for o in objects:
            with self.assertRaises( TypeError ):
                SymmetryOperation( o )

    def test_symmetry_operation_is_initialised_with_label( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        label = 'E'
        so = SymmetryOperation( matrix, label=label ) 
        self.assertEqual( so.label, label )

    def test_mul( self ):
        matrix_a = np.matrix( [ [ 1, 1 ], [ 0, 0 ] ] )
        matrix_b = np.matrix( [ [ 1, 0 ], [ 1, 0 ] ] )
        so = SymmetryOperation( matrix_a )
        np.testing.assert_array_equal( ( so * matrix_b ).matrix , np.matrix( [ [ 2, 0 ], [ 0, 0 ] ] ) )

    def test_mul_with_non_matrix( self ):
        matrix_a = np.matrix( [ [ 1, 1 ], [ 0, 0 ] ] )
        array_b	= np.array( [ [ 1, 0 ], [ 1, 0 ] ] )
        so = SymmetryOperation( matrix_a )
        np.testing.assert_array_equal( ( so * array_b ).matrix, np.matrix( [ [ 2, 0 ], [ 0, 0 ] ] ) )

    def test_invert( self ):
        matrix_a = np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        matrix_b = np.matrix( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] )
        so = SymmetryOperation( matrix_a )
        np.testing.assert_array_equal( so.invert().matrix, matrix_b )
    
    def test_from_vector( self ):
        vector = [ 2, 3, 1 ]
        so = SymmetryOperation.from_vector( vector )
        np.testing.assert_array_equal( so.matrix, np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] ) )    

    def test_from_vector_counting_from_zero( self ):
        vector = [ 1, 2, 0 ]
        so = SymmetryOperation.from_vector( vector, count_from_zero=True )
        np.testing.assert_array_equal( so.matrix, np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] ) )    
if __name__ == '__main__':
    unittest.main()
