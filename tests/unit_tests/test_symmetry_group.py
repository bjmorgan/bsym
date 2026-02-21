import unittest
from bsym import SymmetryGroup, SymmetryOperation
from bsym.configuration import Configuration
from unittest.mock import Mock, patch, call
import numpy as np

class SymmetryGroupTestCase( unittest.TestCase ):
    """Tests for SymmetryGroup class"""

    def test_symmetry_group_is_initialised( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        self.assertEqual( sg.symmetry_operations[0], s0 )
        self.assertEqual( sg.symmetry_operations[1], s1 )

    def test_default_symmetry_operations_are_not_shared(self):
        sg1 = SymmetryGroup()
        sg2 = SymmetryGroup()
        self.assertIsNot(sg1.symmetry_operations, sg2.symmetry_operations)

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

    def test_extend_invalidates_caches(self):
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        # Access to populate caches
        _ = sg.stacked_index_mappings
        _ = sg.unique_index_mappings
        s2 = SymmetryOperation.from_vector([2, 3, 1])
        sg.extend([s2])
        self.assertEqual(sg.stacked_index_mappings.shape[0], 3)

    def test_append( self ):
        s0, s1 = Mock( spec=SymmetryOperation ), Mock( spec=SymmetryOperation )
        sg = SymmetryGroup( symmetry_operations=[ s0, s1 ] )
        s2 = Mock( spec=SymmetryOperation)
        sg.append( s2 )
        self.assertEqual( sg.symmetry_operations, [ s0, s1, s2 ] )

    def test_append_invalidates_caches(self):
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        _ = sg.stacked_index_mappings
        _ = sg.unique_index_mappings
        s2 = SymmetryOperation.from_vector([2, 3, 1])
        sg.append(s2)
        self.assertEqual(sg.stacked_index_mappings.shape[0], 3)
        	     
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

    def test_operate_on(self):
        """Test operate_on with mock operations."""
        s0 = Mock(spec=SymmetryOperation)
        s1 = Mock(spec=SymmetryOperation)
        s0.index_mapping = np.array([0, 1, 2])
        s1.index_mapping = np.array([0, 2, 1])
        configuration = Configuration([1, 0, 0])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        all_configurations = sg.operate_on(configuration)
        self.assertEqual(all_configurations, [Configuration([1, 0, 0]),
                                            Configuration([1, 0, 0])])
        
    def test_operate_on_returns_minimal_set(self):
        """Test operate_on returns minimal set when requested."""
        s0 = Mock(spec=SymmetryOperation)
        s1 = Mock(spec=SymmetryOperation)
        s0.index_mapping = np.array([0, 1, 2])
        s1.index_mapping = np.array([0, 1, 2])  # Same as s0
        configuration = Configuration([1, 0, 0])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        all_configurations = sg.operate_on(configuration, minimal_set=True)
        self.assertEqual(all_configurations, [Configuration([1, 0, 0])])

    def test_stacked_index_mappings_has_correct_shape(self):
        """Test that stacked_index_mappings returns correct shape."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([2, 3, 1])
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        
        mappings = sg.stacked_index_mappings
        self.assertEqual(mappings.shape, (3, 3))  # 3 operations, 3 sites
    
    def test_stacked_index_mappings_contains_correct_values(self):
        """Test that stacked_index_mappings contains correct index mappings."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([2, 3, 1])
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        
        mappings = sg.stacked_index_mappings
        np.testing.assert_array_equal(mappings[0], [0, 1, 2])  # Identity
        np.testing.assert_array_equal(mappings[1], [1, 0, 2])  # Swap first two
        np.testing.assert_array_equal(mappings[2], [2, 0, 1])  # Cycle
    
    def test_stacked_index_mappings_is_cached(self):
        """Test that stacked_index_mappings is cached on repeated access."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        
        mappings1 = sg.stacked_index_mappings
        mappings2 = sg.stacked_index_mappings
        self.assertIs(mappings1, mappings2)  # Same object
    
    def test_unique_index_mappings_removes_duplicates(self):
        """Test that unique_index_mappings removes duplicate operations."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([1, 2, 3])  # Duplicate of s0
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        
        unique_mappings = sg.unique_index_mappings
        self.assertEqual(unique_mappings.shape[0], 2)  # Only 2 unique operations
    
    def test_unique_index_mappings_is_cached(self):
        """Test that unique_index_mappings is cached on repeated access."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        
        mappings1 = sg.unique_index_mappings
        mappings2 = sg.unique_index_mappings
        self.assertIs(mappings1, mappings2)  # Same object
    
    def test_operate_on_batched_returns_same_results_as_original(self):
        """Test that batched operate_on returns same results as original implementation."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([2, 3, 1])
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        config = Configuration([1, 0, 0])
        
        results = sg.operate_on(config, minimal_set=False)
        
        self.assertEqual(len(results), 3)
        np.testing.assert_array_equal(results[0].vector, [1, 0, 0])  # Identity
        np.testing.assert_array_equal(results[1].vector, [0, 1, 0])  # Swap
        np.testing.assert_array_equal(results[2].vector, [0, 1, 0])  # Cycle result
    
    def test_operate_on_minimal_set_removes_duplicates(self):
        """Test that operate_on with minimal_set=True removes duplicates."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([2, 3, 1])
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        
        # Configuration where s0 and s1 give same result
        config = Configuration([1, 1, 0])
        results = sg.operate_on(config, minimal_set=True)
        
        # Should have only 2 unique results
        self.assertEqual(len(results), 2)
        result_tuples = [tuple(r.vector) for r in results]
        self.assertEqual(len(set(result_tuples)), 2)
    
    def test_operate_on_minimal_set_false_includes_all(self):
        """Test that operate_on with minimal_set=False includes all operations."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        
        config = Configuration([1, 1, 0])
        results = sg.operate_on(config, minimal_set=False)
        
        self.assertEqual(len(results), 2)
        
    def test_operate_on_minimal_set_uses_unique_operations(self):
        """Test that operate_on with minimal_set=True only applies unique operations once."""
        # Create a group with 3 operations, where 2 have identical index_mappings
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([1, 2, 3])  # Duplicate of s0
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        
        # Verify we have 3 operations but only 2 unique mappings
        self.assertEqual(len(sg.symmetry_operations), 3)
        self.assertEqual(sg.unique_index_mappings.shape[0], 2)
        
        # When minimal_set=True, should only apply the 2 unique operations
        config = Configuration([1, 0, 0])
        results = sg.operate_on(config, minimal_set=True)
        
        # Should get 2 results (one per unique operation)
        self.assertEqual(len(results), 2)
        
        # When minimal_set=False, should apply all 3 operations
        results_all = sg.operate_on(config, minimal_set=False)
        self.assertEqual(len(results_all), 3)


if __name__ == '__main__':
    unittest.main()
