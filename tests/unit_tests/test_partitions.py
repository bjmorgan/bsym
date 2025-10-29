# tests/unit_tests/test_partitions.py
import unittest
from bsym.partitions import generate_partitions
from bsym.partitions import compute_mapping_vector
from bsym.partitions import satisfies_bounds


class GeneratePartitionsTestCase(unittest.TestCase):
	"""Tests for generate_partitions function"""

	def test_generate_partitions_n4_k2(self):
		"""N=4, k=2 should give 3 partitions"""
		result = list(generate_partitions(4, 2))
		expected = [
			(4, 0),
			(3, 1),
			(2, 2)
		]
		self.assertEqual(result, expected)

	def test_generate_partitions_n4_k3(self):
		"""N=4, k=3 should give 4 partitions"""
		result = list(generate_partitions(4, 3))
		expected = [
			(4, 0, 0),
			(3, 1, 0),
			(2, 2, 0),
			(2, 1, 1)
		]
		self.assertEqual(result, expected)

	def test_generate_partitions_n4_k4(self):
		"""N=4, k=4 should include uniform partition"""
		result = list(generate_partitions(4, 4))
		expected = [
			(4, 0, 0, 0),
			(3, 1, 0, 0),
			(2, 2, 0, 0),
			(2, 1, 1, 0),
			(1, 1, 1, 1)
		]
		self.assertEqual(result, expected)

	def test_partitions_are_tuples(self):
		"""Partitions should be returned as tuples"""
		result = list(generate_partitions(3, 2))
		for partition in result:
			self.assertIsInstance(partition, tuple)

	def test_partitions_have_correct_length(self):
		"""All partitions should have length k"""
		result = list(generate_partitions(5, 3))
		for partition in result:
			self.assertEqual(len(partition), 3)

	def test_partitions_sum_to_n(self):
		"""Each partition should sum to n"""
		n, k = 6, 3
		result = list(generate_partitions(n, k))
		for partition in result:
			self.assertEqual(sum(partition), n)

	def test_partitions_are_descending(self):
		"""Each partition should be sorted descending"""
		result = list(generate_partitions(7, 4))
		for partition in result:
			self.assertEqual(list(partition), sorted(partition, reverse=True))

	def test_generate_partitions_n0_k2(self):
		"""N=0 should give single partition of zeros"""
		result = list(generate_partitions(0, 2))
		expected = [(0, 0)]
		self.assertEqual(result, expected)

	def test_generate_partitions_n1_k1(self):
		"""N=1, k=1 minimal case"""
		result = list(generate_partitions(1, 1))
		expected = [(1,)]
		self.assertEqual(result, expected)

	def test_generate_partitions_n_less_than_k(self):
		"""N<k should work correctly"""
		result = list(generate_partitions(2, 5))
		expected = [
			(2, 0, 0, 0, 0),
			(1, 1, 0, 0, 0)
		]
		self.assertEqual(result, expected)

	def test_no_duplicate_partitions(self):
		"""Should not generate duplicate partitions"""
		result = list(generate_partitions(6, 3))
		self.assertEqual(len(result), len(set(result)))
		
		
class ComputeMappingVectorTestCase(unittest.TestCase):
	"""Tests for compute_mapping_vector function"""

	def test_identity_mapping(self):
		"""Same partition should give identity mapping"""
		result = compute_mapping_vector((2, 1, 1), (2, 1, 1))
		self.assertEqual(result, [0, 1, 2])

	def test_swap_first_two_elements(self):
		"""Swapping first two elements: A₂BC → AB₂C"""
		# Canonical: A=2, B=1, C=1
		# Permuted:  A=1, B=2, C=1
		# Should swap A↔B
		result = compute_mapping_vector((2, 1, 1), (1, 2, 1))
		self.assertEqual(result, [1, 0, 2])

	def test_partitions_with_zeros(self):
		"""Partitions with trailing zeros: A₃B → BA₃"""
		result = compute_mapping_vector((3, 1, 0), (1, 3, 0))
		self.assertEqual(result, [1, 0, 2])

	def test_mapping_is_valid_permutation(self):
		"""Mapping should be a valid permutation of 1..k"""
		result = compute_mapping_vector((3, 2, 1), (1, 3, 2))
		self.assertEqual(sorted(result), list(range(len(result))))
	
		
class SatisfiesBoundsTestCase(unittest.TestCase):
	"""Tests for satisfies_bounds function"""

	def test_composition_within_bounds(self):
		"""Composition satisfying all bounds should pass"""
		composition = {'A': 2, 'B': 1, 'C': 1}
		bounds = {'A': (1, 3), 'B': (1, 2), 'C': (0, 2)}
		self.assertTrue(satisfies_bounds(composition, bounds))

	def test_composition_exceeds_max_bound(self):
		"""Composition exceeding max should fail"""
		composition = {'A': 3, 'B': 1}
		bounds = {'A': (0, 2), 'B': (0, 2)}
		self.assertFalse(satisfies_bounds(composition, bounds))

	def test_composition_below_min_bound(self):
		"""Composition below min should fail"""
		composition = {'A': 1, 'B': 3}
		bounds = {'A': (2, 4), 'B': (0, 4)}
		self.assertFalse(satisfies_bounds(composition, bounds))

	def test_none_bounds_allows_any_value(self):
		"""None in bounds means unbounded"""
		composition = {'A': 10, 'B': 5}
		bounds = {'A': (None, None), 'B': (1, None)}
		self.assertTrue(satisfies_bounds(composition, bounds))

	def test_species_not_in_bounds_defaults_to_allowed(self):
		"""Species not in bounds should be allowed (default 0 to n_sites)"""
		composition = {'A': 2, 'B': 2}
		bounds = {'A': (1, 3)}  # B not specified
		self.assertTrue(satisfies_bounds(composition, bounds))

	def test_species_in_bounds_but_absent_fails_if_min_positive(self):
		"""Species in bounds but not in composition (count=0) must satisfy min"""
		composition = {'A': 4}  # B absent (implicitly 0)
		bounds = {'A': (0, 4), 'B': (1, 4)}  # B requires min=1
		self.assertFalse(satisfies_bounds(composition, bounds))

	def test_species_in_bounds_but_absent_passes_if_min_zero(self):
		"""Absent species OK if min bound is 0 or None"""
		composition = {'A': 4}  # B absent
		bounds = {'A': (0, 4), 'B': (0, 4)}  # B allows 0
		self.assertTrue(satisfies_bounds(composition, bounds))

if __name__ == '__main__':
	unittest.main()