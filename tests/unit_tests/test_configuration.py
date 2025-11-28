import unittest
from unittest.mock import Mock, patch
from bsym.configuration import Configuration
from bsym.configuration import save_configurations, load_configurations
from bsym.symmetry_group import SymmetryGroup
from bsym import SymmetryOperation
import numpy as np
import tempfile
import json
import os

class TestConfiguration(unittest.TestCase):

    def setUp(self):
        self.configuration = Configuration([1, 0, 0])

    def test__eq__returns_true_for_a_match1(self):
        other_configuration = Configuration([1, 0, 0])
        self.assertEqual(self.configuration == other_configuration, True)

    def test__eq__returns_false_for_a_non_match1(self):
        other_configuration = Configuration([1, 1, 0])
        self.assertEqual(self.configuration == other_configuration, False)

    def test_hash_equal_configurations(self):
        config1 = Configuration([1, 1, 0])
        config2 = Configuration([1, 1, 0])
        self.assertEqual(hash(config1), hash(config2))

    def test_hash_different_configurations(self):
        config1 = Configuration([1, 1, 0])
        config2 = Configuration([0, 1, 1])
        self.assertNotEqual(hash(config1), hash(config2))

    def test_hash_consistent(self):
        config = Configuration([1, 1, 0])
        hash_value1 = hash(config)
        hash_value2 = hash(config)
        self.assertEqual(hash_value1, hash_value2)

    def test_matches_true(self):
        config1 = Configuration([1, 1, 0])
        config2 = Configuration([1, 1, 0])
        with patch.object(Configuration, '__eq__', return_value=True):
            self.assertTrue(config1.matches(config2))

    def test_matches_false(self):
        config1 = Configuration([1, 1, 0])
        config2 = Configuration([0, 1, 1])
        with patch.object(Configuration, '__eq__', return_value=False):
            self.assertFalse(config1.matches(config2))

    def test_matches_invalid_type(self):
        config = Configuration([1, 1, 0])
        invalid_config = [1, 1, 0]
        with self.assertRaises(TypeError):
            config.matches(invalid_config)

    def test_is_equivalent_to_if_equivalent(self):
        test_configuration = Configuration([0, 1, 0])
        symmetry_operations = [Mock(spec=SymmetryOperation)]
        symmetry_operations[0].operate_on = Mock(return_value=test_configuration)
        with patch.object(Configuration, '__eq__', return_value=True):
            self.assertEqual(self.configuration.is_equivalent_to(test_configuration, symmetry_operations), True)

    def test_is_equivalent_to_if_not_equivalent( self ):
        test_configuration = Configuration([0, 1, 0 ])
        symmetry_operations = [Mock(spec=SymmetryOperation)] 
        symmetry_operations[0].operate_on = Mock(return_value=Configuration([0, 0, 1]))
        with patch.object(Configuration, '__eq__', return_value=False):
            self.assertEqual(self.configuration.is_equivalent_to(test_configuration, symmetry_operations), False)

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

    def test_map_objects(self):
        self.assertEqual(self.configuration.map_objects(['A', 'B', 'C']), {1: ['A'], 0: ['B', 'C']})  

    def test_map_objects_with_incompatible_object_list_raises_ValueError(self):
        with self.assertRaises(ValueError):
            self.configuration.map_objects(['A', 'B'])
            
    def test_get_byte_equivalents_returns_set_of_bytes(self):
        """Test that get_byte_equivalents returns a set of bytes."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        sg = SymmetryGroup(symmetry_operations=[s0])
        
        config = Configuration([1, 0, 0])
        result = config.get_byte_equivalents(sg)
        
        self.assertIsInstance(result, set)
        for item in result:
            self.assertIsInstance(item, bytes)
    
    def test_get_byte_equivalents_returns_correct_values(self):
        """Test that get_byte_equivalents returns correct byte representations."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        
        config = Configuration([1, 0, 0])
        byte_equivalents = config.get_byte_equivalents(sg)
        
        # Should get byte representations of [1, 0, 0] and [0, 1, 0]
        expected = {
            np.array([1, 0, 0], dtype=np.int8).tobytes(),
            np.array([0, 1, 0], dtype=np.int8).tobytes()
        }
        self.assertEqual(byte_equivalents, expected)
    
    def test_get_byte_equivalents_uses_unique_operations(self):
        """Test that get_byte_equivalents only applies unique operations."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        s1 = SymmetryOperation.from_vector([2, 1, 3])
        s2 = SymmetryOperation.from_vector([1, 2, 3])  # Duplicate of s0
        sg = SymmetryGroup(symmetry_operations=[s0, s1, s2])
        
        config = Configuration([1, 0, 0])
        byte_equivalents = config.get_byte_equivalents(sg)
        
        # Should only apply 2 unique operations, not 3
        self.assertEqual(len(byte_equivalents), 2)
    
    def test_get_byte_equivalents_with_larger_configuration(self):
        """Test get_byte_equivalents with larger configuration space."""
        s0 = SymmetryOperation.from_vector([1, 2, 3, 4, 5])
        s1 = SymmetryOperation.from_vector([5, 4, 3, 2, 1])
        sg = SymmetryGroup(symmetry_operations=[s0, s1])
        
        # Use an asymmetric configuration so reverse gives different result
        config = Configuration([1, 0, 0, 0, 0])
        byte_equivalents = config.get_byte_equivalents(sg)
        
        # Should have 2 results (identity and reverse are different)
        self.assertEqual(len(byte_equivalents), 2)
        # All should be bytes
        for val in byte_equivalents:
            self.assertIsInstance(val, bytes)
    
    def test_get_byte_equivalents_produces_different_values_for_different_configs(self):
        """Test that different configurations produce different byte representations."""
        s0 = SymmetryOperation.from_vector([1, 2, 3])
        sg = SymmetryGroup(symmetry_operations=[s0])
        
        config1 = Configuration([1, 0, 0])
        config2 = Configuration([0, 1, 0])
        
        bytes1 = config1.get_byte_equivalents(sg)
        bytes2 = config2.get_byte_equivalents(sg)
        
        # Different configurations should produce different byte sets
        self.assertNotEqual(bytes1, bytes2)
    
    def test_as_bytes_returns_byte_representation(self):
        """Test that as_bytes returns the byte representation of the configuration."""
        config = Configuration([1, 0, 1])
        result = config.as_bytes()
        
        self.assertIsInstance(result, bytes)
        expected = np.array([1, 0, 1], dtype=np.int8).tobytes()
        self.assertEqual(result, expected)
    
    def test_as_bytes_consistent_with_hash(self):
        """Test that as_bytes is consistent with __hash__."""
        config1 = Configuration([1, 0, 1])
        config2 = Configuration([1, 0, 1])
        
        # Same configuration should have same bytes and hash
        self.assertEqual(config1.as_bytes(), config2.as_bytes())
        self.assertEqual(hash(config1), hash(config2))
        
    def test_tuple_to_bytes_returns_bytes(self):
        """Test that tuple_to_bytes returns bytes."""
        result = Configuration.tuple_to_bytes((1, 0, 1))
        self.assertIsInstance(result, bytes)
    
    def test_tuple_to_bytes_converts_to_int8(self):
        """Test that tuple_to_bytes uses int8 representation."""
        result = Configuration.tuple_to_bytes((1, 0, 1))
        expected = np.array([1, 0, 1], dtype=np.int8).tobytes()
        self.assertEqual(result, expected)
    
    def test_tuple_to_bytes_consistent_results(self):
        """Test that same tuple produces same bytes."""
        result1 = Configuration.tuple_to_bytes((1, 0, 1))
        result2 = Configuration.tuple_to_bytes((1, 0, 1))
        self.assertEqual(result1, result2)
    
    def test_tuple_to_bytes_different_for_different_tuples(self):
        """Test that different tuples produce different bytes."""
        result1 = Configuration.tuple_to_bytes((1, 0, 0))
        result2 = Configuration.tuple_to_bytes((0, 1, 0))
        self.assertNotEqual(result1, result2)
    
    def test_array_to_bytes_returns_bytes(self):
        """Test that array_to_bytes returns bytes."""
        arr = np.array([1, 0, 1], dtype=np.int8)
        result = Configuration.array_to_bytes(arr)
        self.assertIsInstance(result, bytes)
    
    def test_array_to_bytes_uses_array_directly(self):
        """Test that array_to_bytes converts int8 array to bytes."""
        arr = np.array([1, 0, 1], dtype=np.int8)
        result = Configuration.array_to_bytes(arr)
        expected = arr.tobytes()
        self.assertEqual(result, expected)
    
    def test_array_to_bytes_consistent_with_tuple_to_bytes(self):
        """Test that tuple and array conversions produce same bytes for same values."""
        tup = (1, 0, 1)
        arr = np.array([1, 0, 1], dtype=np.int8)
        
        tuple_result = Configuration.tuple_to_bytes(tup)
        array_result = Configuration.array_to_bytes(arr)
        
        self.assertEqual(tuple_result, array_result)
    
    def test_as_bytes_uses_array_to_bytes(self):
        """Test that as_bytes is consistent with array_to_bytes."""
        config = Configuration([1, 0, 1])
        
        instance_result = config.as_bytes()
        static_result = Configuration.array_to_bytes(config.vector)
        
        self.assertEqual(instance_result, static_result)
        
    def test_configuration_to_dict(self):
        """
        Test that to_dict returns a JSON-serialisable dictionary.
        """
        config = Configuration([0, 1, 0, 1])
        result = config.to_dict()
        self.assertEqual(result, {'vector': [0, 1, 0, 1]})
    
    def test_configuration_from_dict(self):
        """
        Test that from_dict creates a Configuration from a dictionary.
        """
        d = {'vector': [1, 0, 1, 0]}    
        config = Configuration.from_dict(d)
        self.assertEqual(config.tolist(), [1, 0, 1, 0])
    
    def test_configuration_round_trip(self):
        """
        Test that to_dict and from_dict round-trip correctly.
        """
        original = Configuration([0, 0, 1, 1])    
        d = original.to_dict()
        restored = Configuration.from_dict(d)
        self.assertEqual(original.tolist(), restored.tolist())
        
    def test_save_configurations_creates_json_file(self):
        """
        Test that save_configurations creates a JSON file.
        """
        configs = [Configuration([0, 1, 0, 1]), Configuration([1, 0, 1, 0])]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            filename = f.name
        
        try:
            save_configurations(configs, filename)
            
            with open(filename) as f:
                data = json.load(f)
            
            self.assertEqual(len(data), 2)
            self.assertEqual(data[0], {'vector': [0, 1, 0, 1]})
            self.assertEqual(data[1], {'vector': [1, 0, 1, 0]})
        finally:
            os.unlink(filename)

    def test_load_configurations_reads_json_file(self):
        """
        Test that load_configurations reads configurations from a JSON file.
        """
        data = [{'vector': [0, 1, 0, 1]}, {'vector': [1, 0, 1, 0]}]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(data, f)
            filename = f.name
        
        try:
            configs = load_configurations(filename)
            
            self.assertEqual(len(configs), 2)
            self.assertEqual(configs[0].tolist(), [0, 1, 0, 1])
            self.assertEqual(configs[1].tolist(), [1, 0, 1, 0])
        finally:
            os.unlink(filename)
    
    def test_save_and_load_configurations_round_trip(self):
        """
        Test that save and load round-trip correctly.
        """
        original_configs = [
            Configuration([0, 0, 1, 1]),
            Configuration([0, 1, 0, 1]),
            Configuration([1, 1, 0, 0]),
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            filename = f.name
        
        try:
            save_configurations(original_configs, filename)
            loaded_configs = load_configurations(filename)
            
            self.assertEqual(len(loaded_configs), len(original_configs))
            for orig, loaded in zip(original_configs, loaded_configs):
                self.assertEqual(orig.tolist(), loaded.tolist())
        finally:
            os.unlink(filename)

if __name__ == '__main__':
    unittest.main()
