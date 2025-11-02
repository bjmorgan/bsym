# Composition Enumeration

## The Composition Enumeration Problem

### Single vs Multiple Compositions

The `unique_configurations()` method finds all symmetry-inequivalent configurations for a **fixed composition** - a specific distribution of object types. For example, `{0: 2, 1: 2}` specifies exactly 2 objects of type 0 and 2 of type 1.

However, many problems require exploring **multiple compositions**:

- Building phase diagrams (e.g., Li<sub>x</sub>Na<sub>1-x</sub>CoO<sub>2</sub> across all x)
- Finding energetically favorable compositions
- Surveying disorder across different stoichiometries

### The Challenge

For an n-site system with k species, there are many possible compositions. A naive approach would call `unique_configurations()` separately for each composition, but this is inefficient because many compositions are related by symmetry.

### Example: Binary System

Consider a 4-site system with 2 species (labels 0 and 1):

Possible compositions:
- `(4, 0)` - all sites are species 0
- `(3, 1)` - 3 sites species 0, 1 site species 1
- `(2, 2)` - 2 sites of each species
- `(1, 3)` - 1 site species 0, 3 sites species 1
- `(0, 4)` - all sites are species 1

If we are studying Li/Na substitution on 4 sites, `(3, 1)` might represent Li<sub>3</sub>Na, while `(1, 3)` represents LiNa<sub>3</sub>. These are different chemically, but from a **configuration enumeration perspective** they are equivalent &mdash; the same counting problem with labels swapped.

## Species Exchange Symmetry

### What is Species Exchange Symmetry?

Two compositions are related by **species exchange symmetry** if one can be obtained from the other by permuting the species labels.

For a composition represented as a tuple `(n₀, n₁, n₂, ...)` where n<sub>i</sub> is the count of species i:

- `(3, 1, 0)` and `(1, 3, 0)` are related by swapping species 0 ↔ species 1
- `(2, 1, 1)` and `(1, 2, 1)` are related by swapping species 0 ↔ species 1
- `(2, 1, 1)`, `(1, 2, 1)`, `(1, 1, 2)` are all related by various species permutations

### Why This Matters

Compositions related by species exchange produce **identical sets** of symmetry-inequivalent configurations (up to relabeling).

If we find the unique configurations for composition `(3, 1)`, we can generate the configurations for `(1, 3)` simply by swapping the species labels in each configuration. No expensive symmetry analysis is needed for the second composition.

## Partitions and Canonical Compositions

### Integer Partitions

An **integer partition** of n into k parts is a way of writing n as a sum of k non-negative integers, in non-increasing order.

For n=4 sites and k=2 species:
- `(4, 0)` - partition: 4 + 0 = 4
- `(3, 1)` - partition: 3 + 1 = 4
- `(2, 2)` - partition: 2 + 2 = 4

Note: `(3, 1)` and `(1, 3)` represent the **same partition**. When written in non-increasing order, both become `(3, 1)`.

### Partitions Group Compositions

Each partition corresponds to all compositions related by species exchange:

**Partition (3, 1)** corresponds to:
- Composition `(3, 1)` - 3 of species 0, 1 of species 1
- Composition `(1, 3)` - 1 of species 0, 3 of species 1

**Partition (2, 2)** corresponds to:
- Composition `(2, 2)` - 2 of each species (only one arrangement - symmetric)

**Partition (4, 0)** corresponds to:
- Composition `(4, 0)` - all sites species 0
- Composition `(0, 4)` - all sites species 1

### Canonical Composition

The **canonical composition** for a partition is the partition tuple itself.

- For partition `(3, 1)`: canonical composition is `(3, 1)`
- For partition `(2, 2)`: canonical composition is `(2, 2)`

All other permutations of the partition are **non-canonical compositions**.

### The Optimization Strategy

Instead of performing symmetry analysis on every composition:

1. Generate all partitions of n sites into k species
2. For each partition, perform symmetry analysis **only on the canonical composition**
3. For non-canonical compositions, generate configurations by **relabeling species** in the canonical results

This reduces computational cost from O(compositions) to O(partitions).

### Example: 4 Sites, 2 Species

All compositions: `(4,0)`, `(3,1)`, `(2,2)`, `(1,3)`, `(0,4)` - **5 total**

All partitions: `(4,0)`, `(3,1)`, `(2,2)` - **3 total**

**Without optimization** (naive approach):
- 5 separate symmetry analyses (one per composition)

**With optimization**:
- Analyze partition `(4,0)` → generate composition `(0,4)` by relabeling
- Analyze partition `(3,1)` → generate composition `(1,3)` by relabeling
- Analyze partition `(2,2)` (symmetric - no additional compositions)
- **3 symmetry analyses**

**Reduction: 40%** (5 analyses → 3 analyses)

## Species Mapping and Relabeling

### The Mapping Vector

To generate configurations for a non-canonical composition from canonical results, we compute a **species mapping vector** that describes how to relabel species.

For partition `(3, 1)` with compositions `(3, 1)` → `(1, 3)`:
- Species 0 in canonical (count 3) appears at position 1 in `(1, 3)` → mapping[0] = 1
- Species 1 in canonical (count 1) appears at position 0 in `(1, 3)` → mapping[1] = 0

Mapping vector: `[1, 0]` (swap species 0↔1)

### Applying the Mapping

Apply the mapping vector to each configuration by relabeling species:

**Canonical configuration for (3, 1)**:
```
[0, 0, 0, 1]  # 3 sites are species 0, 1 site is species 1
```

**After applying mapping [1, 0]**:
```
[1, 1, 1, 0]  # 3 sites are species 1, 1 site is species 0
```

This produces a valid configuration for composition `(1, 3)`.

The degeneracy (count) is preserved during relabeling - if the canonical configuration represents 4 equivalent arrangements, the relabeled configuration also represents 4 equivalent arrangements.

## The Algorithm

### Overview

The composition enumeration algorithm:

1. **Generate partitions**: Find all partitions of n sites into k species
2. **Generate compositions**: For each partition, create all permutations to get all compositions
3. **Filter by bounds** (if specified): Keep only compositions satisfying occupancy constraints
4. **Analyze canonical compositions**: For each partition, perform symmetry analysis on the canonical composition only
5. **Generate non-canonical results**: For other compositions from the same partition, relabel the canonical configurations using species mapping vectors

### Step-by-Step Process

**Given**: n-site configuration space, k species, optional occupancy bounds

**Step 1: Generate Partitions**

Generate all integer partitions of n into k parts:
```python
# For n=4, k=2:
partitions = [(4, 0), (3, 1), (2, 2)]
```

**Step 2: Generate Compositions from Each Partition**

For each partition, generate all distinct permutations:
```python
# For partition (3, 1):
permutations = [(3, 1), (1, 3)]
```

**Step 3: Filter by Bounds (Optional)**

If occupancy bounds are specified, check each composition:
```python
# bounds = {0: (1, 3), 1: (1, 3)}
# Composition (4, 0) violates bounds (species 1 has 0, needs at least 1)
# Composition (3, 1) satisfies bounds ✓
```

**Step 4: Analyze Canonical Composition**

For each partition, perform expensive symmetry analysis on the canonical composition:
```python
# For partition (3, 1):
canonical = (3, 1)
site_distribution = {0: 3, 1: 1}
canonical_configs = config_space.unique_configurations(site_distribution)
```

**Step 5: Generate Non-Canonical Configurations**

For each non-canonical composition from the same partition:
```python
# For composition (1, 3) from partition (3, 1):
mapping = compute_mapping_vector((3, 1), (1, 3))  # Returns [1, 0]
relabeled_configs = [apply_species_mapping(config, mapping) 
                     for config in canonical_configs]
results[(1, 3)] = relabeled_configs
```

### Implementation in bsym

The algorithm is implemented in `ConfigurationSpace.unique_configurations_by_composition()`:
```python
from bsym import ConfigurationSpace

config_space = ConfigurationSpace(
    objects=[1, 2, 3, 4],
    symmetry_group=my_symmetry_group
)

# Enumerate all compositions for 2 species
results = config_space.unique_configurations_by_composition(n_species=2)

# Results is a dictionary: {composition_tuple: [Configuration, ...]}
# e.g., {(4, 0): [...], (3, 1): [...], (2, 2): [...], (1, 3): [...], (0, 4): [...]}
```

With optional bounds:
```python
# Only compositions with 1-3 of each species
bounds = {0: (1, 3), 1: (1, 3)}
results = config_space.unique_configurations_by_composition(
    n_species=2,
    bounds=bounds
)

# Results: {(3, 1): [...], (2, 2): [...], (1, 3): [...]}
# (4, 0) and (0, 4) excluded by bounds
```
