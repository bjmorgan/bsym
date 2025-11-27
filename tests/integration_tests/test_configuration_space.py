import unittest
from bsym import ConfigurationSpace, SymmetryGroup, SymmetryOperation


class TestRandomUniqueConfigurations(unittest.TestCase):
	"""Integration tests for ConfigurationSpace.random_unique_configurations."""

	def setUp(self):
		"""Set up C4v symmetry group and configuration space."""
		e = SymmetryOperation.from_vector([1, 2, 3, 4], label='E')
		c4 = SymmetryOperation.from_vector([2, 3, 4, 1], label='C4')
		c4_inv = SymmetryOperation.from_vector([4, 1, 2, 3], label='C4i')
		c2 = SymmetryOperation.from_vector([3, 4, 1, 2], label='C2')
		sigma_x = SymmetryOperation.from_vector([4, 3, 2, 1], label='s_x')
		sigma_y = SymmetryOperation.from_vector([2, 1, 4, 3], label='s_y')
		sigma_ac = SymmetryOperation.from_vector([1, 4, 3, 2], label='s_ac')
		sigma_bd = SymmetryOperation.from_vector([3, 2, 1, 4], label='s_bd')
		
		self.c4v = SymmetryGroup([e, c4, c4_inv, c2, sigma_x, sigma_y, sigma_ac, sigma_bd])
		
		self.config_space = ConfigurationSpace(
			objects=['a', 'b', 'c', 'd'],
			symmetry_group=self.c4v
		)

	def test_returned_configurations_are_mutually_inequivalent(self):
		"""
		Test that all returned configurations are symmetry-inequivalent
		to each other.
		"""
		result = self.config_space.random_unique_configurations(
			site_distribution={1: 2, 0: 2},
			n=2,
			seed=42,
		)
		
		for i, config_i in enumerate(result):
			for j, config_j in enumerate(result):
				if i >= j:
					continue
				self.assertFalse(
					config_i.is_equivalent_to(config_j, self.c4v.symmetry_operations),
					f"Configuration {i} and {j} are equivalent"
				)

	def test_same_seed_produces_same_results(self):
		"""
		Test that using the same seed produces identical results.
		"""
		result_1 = self.config_space.random_unique_configurations(
			site_distribution={1: 2, 0: 2},
			n=2,
			seed=42,
		)
		
		result_2 = self.config_space.random_unique_configurations(
			site_distribution={1: 2, 0: 2},
			n=2,
			seed=42,
		)
		
		self.assertEqual(len(result_1), len(result_2))
		for config_1, config_2 in zip(result_1, result_2):
			self.assertEqual(config_1.tolist(), config_2.tolist())
			
	def test_degeneracies_are_correct(self):
		"""
		Test that returned configurations have correct degeneracy values.
		"""
		# For C4v with {1: 2, 0: 2}, there are exactly 2 unique configurations:
		# - adjacent sites: degeneracy 4
		# - diagonal sites: degeneracy 2
		result = self.config_space.random_unique_configurations(
			site_distribution={1: 2, 0: 2},
			n=2,
			seed=42,
		)
		
		degeneracies = sorted([config.count for config in result])
		self.assertEqual(degeneracies, [2, 4])
		
	def test_uniform_sampling_distribution_differs_from_degeneracy_weighted(self):
		"""
		Test that uniform sampling produces a different distribution
		than degeneracy_weighted sampling.
		
		With degeneracy_weighted, high-degeneracy configurations are more
		likely to be sampled. With uniform, all equivalence classes have
		equal probability.
		"""
		n_samples = 100
		
		# Count how often we get the high-degeneracy (adjacent) configuration
		# vs low-degeneracy (diagonal) configuration
		
		degeneracy_weighted_high_count = 0
		for i in range(n_samples):
			result = self.config_space.random_unique_configurations(
				site_distribution={1: 2, 0: 2},
				n=1,
				sampling='degeneracy_weighted',
				seed=i,
			)
			if result[0].count == 4:  # High degeneracy (adjacent)
				degeneracy_weighted_high_count += 1
		
		uniform_high_count = 0
		for i in range(n_samples):
			result = self.config_space.random_unique_configurations(
				site_distribution={1: 2, 0: 2},
				n=1,
				sampling='uniform',
				seed=i,
			)
			if result[0].count == 4:  # High degeneracy (adjacent)
				uniform_high_count += 1
		
		# With degeneracy_weighted: P(adjacent) = 4/6 ≈ 0.67
		# With uniform: P(adjacent) = 0.5
		# So degeneracy_weighted should have more high-degeneracy samples
		
		self.assertGreater(
			degeneracy_weighted_high_count,
			uniform_high_count,
			"degeneracy_weighted should favour high-degeneracy configurations"
		)