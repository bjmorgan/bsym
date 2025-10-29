from typing import TypeVar

K = TypeVar('K', int, str)

def generate_partitions(
	n: int,
	k: int,
	max_value: int | None = None
) -> list[tuple[int, ...]]:
	"""
	Generate all partitions of n into at most k parts.
	
	Args:
		n: Number to partition
		k: Maximum number of parts
		max_value: Maximum value for any part (for recursion)
	"""
	if max_value is None:
		max_value = n  # First call: no restriction
	
	# Base cases
	if k == 0:
		return [()] if n == 0 else [] # can only partition 0 into 0 parts
	if n == 0:
		return [(0,) * k] # pad with zeros
	
	result = []
	
	# Try each possible first part from max_value down to 0
	for first in range(min(n, max_value), -1, -1):
		# Recursively partition remainder
		for rest in generate_partitions(n - first, # remainder to partition
										k - 1,     # remaining number of parts,
										first):    # ensures parts in rest are ≤ first
			result.append((first,) + rest)
	
	return result
	
def compute_mapping_vector(
	canonical_partition: tuple[int, ...],
	permuted_partition: tuple[int, ...]
) -> list[int]:
	"""
	Compute 0-indexed mapping vector from canonical to permuted partition.
	
	Args:
		canonical_partition: The canonical (first) permutation of a partition.
		permuted_partition: A permutation of the canonical partition.
		
	Returns:
		A 0-indexed list where mapping[i] indicates which species position
		in the permuted partition corresponds to position i in canonical.
		
	Example:
		>>> compute_mapping_vector((2, 1, 1), (1, 2, 1))
		[1, 0, 2]  # Swap species 0↔1, keep species 2
	"""
	mapping = []
	permuted_list = list(permuted_partition)
	used = [False] * len(permuted_list)
	
	for canonical_val in canonical_partition:
		# Find matching value in permuted that hasn't been used
		for i, permuted_val in enumerate(permuted_list):
			if permuted_val == canonical_val and not used[i]:
				mapping.append(i)  # Changed from i + 1
				used[i] = True
				break
	
	return mapping

def satisfies_bounds(
	composition: dict[K, int],
	occupancy_bounds: dict[K, tuple[int | None, int | None]]
) -> bool:
	"""
	Check if composition satisfies occupancy bounds.
	
	Args:
		composition: e.g., {'A': 2, 'B': 1, 'C': 1}
		occupancy_bounds: e.g., {'A': (1, 3), 'B': (0, 2)}
		
	Notes:
		Species not in occupancy_bounds default to (0, n_sites).
		Species in occupancy_bounds but not in composition are treated as count=0.
	"""
	n_sites = sum(composition.values())
	
	# Check all species mentioned in either composition or bounds
	all_species = set(composition.keys()) | set(occupancy_bounds.keys())
	
	for species in all_species:
		count = composition.get(species, 0)  # 0 if absent from composition
		lower, upper = occupancy_bounds.get(species, (0, n_sites))
		lower = lower if lower is not None else 0
		upper = upper if upper is not None else n_sites
		if count < lower or count > upper:
			return False
	return True