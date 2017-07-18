import unittest
import numpy as np
from bsym import SymmetryOperation
from bsym.configuration import Configuration
from unittest.mock import patch
import io

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

    def test_from_vector_with_label( self ):
        vector = [ 2, 3, 1 ]
        label = 'A'
        so = SymmetryOperation.from_vector( vector, label=label )
        np.testing.assert_array_equal( so.matrix, np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] ) )
        self.assertEqual( so.label, label )

    def test_symmetry_operation_is_initialised_with_label( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        label = 'E'
        so = SymmetryOperation( matrix, label=label )
        self.assertEqual( so.label, label )

    def test_from_vector_counting_from_zero( self ):
        vector = [ 1, 2, 0 ]
        so = SymmetryOperation.from_vector( vector, count_from_zero=True )
        np.testing.assert_array_equal( so.matrix, np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] ) )    

    def test_similarity_transform( self ):
        matrix_a = np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        matrix_b = np.matrix( [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ] )
        matrix_c = np.linalg.inv( matrix_a )
        so_a = SymmetryOperation( matrix_a )
        so_b = SymmetryOperation( matrix_b )
        np.testing.assert_array_equal( so_a.similarity_transform( so_b ).matrix, matrix_c )

    def test_operate_on( self ):
        matrix = np.matrix( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        so = SymmetryOperation( matrix )
        configuration = np.matrix( [ [ 1, 1, 0 ] ] ).T
        with patch( 'bsym.configuration.Configuration' ) as mock_configuration:
            mock_configuration.return_value = 'foo'
            so.operate_on( configuration )
            self.assertEqual( mock_configuration.call_args[0][0], ( so.matrix * configuration ).tolist() )
    def test_character( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        self.assertEqual( so.character(), 2 )

    def test_as_vector( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        self.assertEqual( so.as_vector(), [ 1, 2 ] )
  
    def test_as_vector_counting_from_zero( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        self.assertEqual( so.as_vector( count_from_zero=True ), [ 0, 1 ] )

    def test_se_label( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        so.set_label( 'new_label' )
        self.assertEqual( so.label, 'new_label' )

    def test_pprint( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        with patch( 'sys.stdout', new=io.StringIO() ) as mock_stdout:
            so.pprint()
            self.assertEqual( mock_stdout.getvalue(), '--- : 1 2\n' ) 

    def test_pprint_with_label( self ):
        matrix = np.matrix( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix, label='L' )
        with patch( 'sys.stdout', new=io.StringIO() ) as mock_stdout:
            so.pprint()
            self.assertEqual( mock_stdout.getvalue(), 'L : 1 2\n' ) 

if __name__ == '__main__':
    unittest.main()
