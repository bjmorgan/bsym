import unittest
from bsym import SymmetryGroup, SymmetryOperation
from unittest.mock import Mock, patch, call
import numpy as np

class SymmetryGroupTestCase( unittest.TestCase ):
    """Tests for SymmetryGroup class"""

    def test_symmetry_group_is_initialised( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        self.assertEqual( sg.symmetry_operations[0], s0 )
        self.assertEqual( sg.symmetry_operations[1], s1 )

    def test_read_from_file( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        with patch( 'numpy.loadtxt' ) as mock_np_loadtxt:
            mock_np_loadtxt.return_value = np.array( [ [ 1, 2 ], [ 2, 1 ] ] )
            with patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' ) as mock_from_vector:
                mock_from_vector.side_effect = [ s0, s1 ]
                sg = SymmetryGroup.read_from_file( 'mock_filename' )
                self.assertEqual( sg.symmetry_operations[0], s0 )
                self.assertEqual( sg.symmetry_operations[1], s1 )
                self.assertEqual( mock_from_vector.call_args_list[0], call( [ 1, 2 ] ) )
                self.assertEqual( mock_from_vector.call_args_list[1], call( [ 2, 1 ] ) )

    def test_read_from_file_with_labels( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        with patch( 'numpy.genfromtxt' ) as mock_np_genfromtxt:
            mock_np_genfromtxt.return_value = np.array( [ [ 'E', '1', '2' ], [ 'C2', '2', '1' ] ] )
            with patch( 'bsym.symmetry_operation.SymmetryOperation.from_vector' ) as mock_from_vector:
                mock_from_vector.side_effect = [ s0, s1 ]
                sg = SymmetryGroup.read_from_file_with_labels( 'mock_filename' )
                self.assertEqual( sg.symmetry_operations[0], s0 )
                self.assertEqual( sg.symmetry_operations[1], s1 )
                self.assertEqual( mock_from_vector.call_args_list[0], call( [ 1, 2 ] ) )
                self.assertEqual( mock_from_vector.call_args_list[1], call( [ 2, 1 ] ) )
                self.assertEqual( s0.set_label.call_args, call( 'E' ) )
                self.assertEqual( s1.set_label.call_args, call( 'C2' ) )
    
    def test_save_symmetry_operation_vectors_to( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        s0.as_vector.return_value = [ 1, 2 ]
        s1.as_vector.return_value = [ 2, 1 ]
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        with patch( 'numpy.savetxt' ) as mock_savetxt:
            sg.save_symmetry_operation_vectors_to( 'filename' ) 
            self.assertEqual( mock_savetxt.call_args[0][0], 'filename' )
            np.testing.assert_array_equal( mock_savetxt.call_args[0][1], np.array( [ [ 1, 2 ], [ 2, 1 ] ] ) )
    
    def test_extend( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        s2 = Mock( spec=SymmetryOperation)
        sg.extend( [ s2 ] )
        self.assertEqual( sg.symmetry_operations, [ s0, s1, s2 ] )
       
    def test_append( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        s2 = Mock( spec=SymmetryOperation)
        sg.append( s2 )
        self.assertEqual( sg.symmetry_operations, [ s0, s1, s2 ] )
        	     
    def test_by_label( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        s0.label = 'A'
        s1.label = 'B'
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        self.assertEqual( sg.by_label( 'A' ), s0 )
        self.assertEqual( sg.by_label( 'B' ), s1 )
  
    def test_labels( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        s0.label = 'A'
        s1.label = 'B'
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        self.assertEqual( sg.labels, [ 'A', 'B' ] )
  
if __name__ == '__main__':
    unittest.main()
