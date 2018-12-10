import unittest
import numpy as np
from bsym import SymmetryOperation, Configuration
from unittest.mock import patch
import io
from bsym.symmetry_operation import is_square, is_permutation_matrix

class SymmetryOperationTestCase( unittest.TestCase ):
    """Tests for symmetry operation functions"""

    def test_symmetry_operation_is_initialised_from_a_matrix( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        np.testing.assert_array_equal( so.matrix, matrix )

    def test_symmetry_operation_is_initialised_from_an_array( self ):
        array = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( array )
        np.testing.assert_array_equal( so.matrix, np.array( array ) )

    def test_symmetry_operation_is_initialised_from_a_list( self ):
        this_list = [ [ 1, 0 ], [ 0, 1 ] ]
        so = SymmetryOperation( this_list )
        np.testing.assert_array_equal( so.matrix, np.array( this_list ) )

    def test_symmetry_operation_raises_typeerror_for_invalid_type( self ):
        objects = [ 'foo', 1, None ]
        for o in objects:
            with self.assertRaises( TypeError ):
                SymmetryOperation( o )

    def test_symmetry_operation_raises_valueerror_for_nonsquare_matrix( self ):
        array = np.array( [ [ 1, 0, 0 ], [ 0, 0, 1 ] ] )
        with self.assertRaises( ValueError ):
            SymmetryOperation( array )

    def test_symmetry_operation_is_initialised_with_label( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        label = 'E'
        so = SymmetryOperation( matrix, label=label ) 
        self.assertEqual( so.label, label )

    def test_mul( self ):
        matrix_a = np.array( [ [ 0, 1 ], [ 1, 0 ] ] )
        matrix_b = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so_a = SymmetryOperation( matrix_a )
        so_b = SymmetryOperation( matrix_b )
        np.testing.assert_array_equal( ( so_a * so_b ).matrix , np.array( [ [ 0, 1 ], [ 1, 0 ] ] ) )

    def test_mul_with_configuration( self ):
        so = SymmetryOperation.from_vector( [ 2, 3, 1 ] )
        conf = Configuration( [ 1, 2, 3 ] )
        new_conf = so * conf
        self.assertEqual( type( new_conf ), Configuration )
        self.assertEqual( new_conf.matches( Configuration( [ 3, 1, 2 ] ) ), True )

    def test_mul_raises_TypeError_with_invalid_type( self ):
        so = SymmetryOperation.from_vector( [ 2, 3, 1 ] )
        with self.assertRaises( TypeError ):
            new_conf = so * 'foo'

    def test_invert( self ):
        matrix_a = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        matrix_b = np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] )
        so = SymmetryOperation( matrix_a )
        np.testing.assert_array_equal( so.invert().matrix, matrix_b )

    def test_invert_sets_label( self ):
        matrix_a = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        so = SymmetryOperation( matrix_a ).invert( label='A' )
        self.assertEqual( so.label, 'A' )
    
    def test_from_vector( self ):
        vector = [ 2, 3, 1 ]
        so = SymmetryOperation.from_vector( vector )
        np.testing.assert_array_equal( so.matrix, np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )    

    def test_from_vector_with_label( self ):
        vector = [ 2, 3, 1 ]
        label = 'A'
        so = SymmetryOperation.from_vector( vector, label=label )
        np.testing.assert_array_equal( so.matrix, np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )
        self.assertEqual( so.label, label )

    def test_symmetry_operation_is_initialised_with_label( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        label = 'E'
        so = SymmetryOperation( matrix, label=label )
        self.assertEqual( so.label, label )

    def test_from_vector_counting_from_zero( self ):
        vector = [ 1, 2, 0 ]
        so = SymmetryOperation.from_vector( vector, count_from_zero=True )
        np.testing.assert_array_equal( so.matrix, np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] ) )    

    def test_similarity_transform( self ):
        matrix_a = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        matrix_b = np.array( [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ] )
        matrix_c = np.linalg.inv( matrix_a )
        so_a = SymmetryOperation( matrix_a )
        so_b = SymmetryOperation( matrix_b )
        np.testing.assert_array_equal( so_a.similarity_transform( so_b ).matrix, matrix_c )

    def test_similarity_transform_with_label( self ):
        matrix_a = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        matrix_b = np.array( [ [ 1, 0, 0 ], [ 0, 0, 1 ], [ 0, 1, 0 ] ] )
        matrix_c = np.linalg.inv( matrix_a )
        so_a = SymmetryOperation( matrix_a )
        so_b = SymmetryOperation( matrix_b )
        label = 'foo'
        np.testing.assert_array_equal( so_a.similarity_transform( so_b, label=label ).label, label )

    def test_operate_on( self ):
        matrix = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        so = SymmetryOperation( matrix )
        configuration = Configuration( [ 1, 2, 3 ] )
        so.operate_on( configuration )
        np.testing.assert_array_equal( so.operate_on( configuration ).vector, np.array( [ 2, 3, 1 ] ) )  

    def test_operate_on_raises_TypeError_with_invalid_type( self ):
        matrix = np.array( [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] )
        so = SymmetryOperation( matrix )
        with self.assertRaises( TypeError ):
            so.operate_on( 'foo' )

    def test_character( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        self.assertEqual( so.character(), 2 )

    def test_as_vector( self ):
        matrix = np.array( [ [ 0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ] )
        so = SymmetryOperation( matrix )
        self.assertEqual( so.as_vector(), [ 2, 3, 1 ] )
  
    def test_as_vector_counting_from_zero( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        self.assertEqual( so.as_vector( count_from_zero=True ), [ 0, 1 ] )

    def test_se_label( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        so.set_label( 'new_label' )
        self.assertEqual( so.label, 'new_label' )

    def test_pprint( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix )
        with patch( 'sys.stdout', new=io.StringIO() ) as mock_stdout:
            so.pprint()
            self.assertEqual( mock_stdout.getvalue(), '--- : 1 2\n' ) 

    def test_pprint_with_label( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix, label='L' )
        with patch( 'sys.stdout', new=io.StringIO() ) as mock_stdout:
            so.pprint()
            self.assertEqual( mock_stdout.getvalue(), 'L : 1 2\n' ) 

    def test_repr( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        so = SymmetryOperation( matrix, label='L' )
        this_repr = so.__repr__()
        self.assertNotEqual( this_repr.find( 'L' ), 0 )
        self.assertNotEqual( this_repr.find( "[[1, 0],\n[0, 1]]" ), 0 )

class SymmetryOperationModuleFunctionsTestCase( unittest.TestCase ):

    def test_is_square_returns_true_if_matrix_is_square( self ):
        matrix = np.array( [ [ 1, 0 ], [ 0, 1 ] ] )
        self.assertEqual( is_square( matrix ), True )

    def test_is_square_returns_false_if_matrix_is_not_square( self ):
        matrix = np.array( [ [ 1, 0, 1 ], [ 0, 1, 1 ] ] )
        self.assertEqual( is_square( matrix ), False )

    def test_is_permutation_matrix_returns_true_if_true( self ):
        matrix = np.array( [ [ 0, 1 ], [ 1, 0 ] ] )
        self.assertEqual( is_permutation_matrix( matrix ), True ) 

    def test_is_permutation_matrix_returns_false_if_false( self ):
        matrix = np.array( [ [ 1, 1 ], [ 0, 0 ] ] )
        self.assertEqual( is_permutation_matrix( matrix ), False ) 



    
if __name__ == '__main__':
    unittest.main()
