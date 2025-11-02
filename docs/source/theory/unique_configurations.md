# Unique Configuration Enumeration

## Configuration Equivalence

Two configurations are **equivalent** if one can be transformed into the other by applying a symmetry operation from the configuration space's symmetry group.

### Example: Square Configuration Space

Consider four sites arranged in a square:

![Square configuration space](./figures/square_configuration_space.pdf)

Without considering symmetry, there are six distinct ways to place two occupied sites (label 1) and two unoccupied sites (label 0):

- `[0, 0, 1, 1]` - positions 2 and 3 occupied
- `[0, 1, 0, 1]` - positions 1 and 3 occupied  
- `[0, 1, 1, 0]` - positions 1 and 2 occupied
- `[1, 0, 0, 1]` - positions 0 and 3 occupied
- `[1, 0, 1, 0]` - positions 0 and 2 occupied
- `[1, 1, 0, 0]` - positions 0 and 1 occupied

However, if the square has C<sub>4v</sub> symmetry (rotations and reflections), some of these configurations become equivalent. For instance, `[0, 0, 1, 1]` (adjacent sites occupied) is equivalent to `[0, 1, 1, 0]`, `[1, 1, 0, 0]`, and `[1, 0, 0, 1]` through rotations.

### Checking Equivalence

To check if configuration A is equivalent to configuration B under a symmetry group:

1. Apply each symmetry operation in the group to configuration A
2. If any resulting configuration equals B, then A and B are equivalent
3. If no operation transforms A into B, they are inequivalent

In code:
```python
from bsym import Configuration, ConfigurationSpace, SymmetryOperation

# Create two configurations
config_a = Configuration([0, 0, 1, 1])
config_b = Configuration([1, 1, 0, 0])

# Check if equivalent under a rotation
c2 = SymmetryOperation.from_vector([3, 4, 1, 2], label='C2')
result = c2.operate_on(config_a)

print(result == config_b)  # True - they are equivalent
```

## The Enumeration Problem

### Finding All Unique Configurations

Given:
- A configuration space with n positions
- A symmetry group with symmetry operations
- A site distribution specifying how many objects of each type to place

The goal is to find **all symmetry-inequivalent configurations**. This is the complete set of configurations where no two members are equivalent to each other under the symmetry group.

Equivalently, this is the minimal set that represents all possible arrangements: every possible permutation of the site distribution is equivalent (under symmetry) to exactly one configuration in this set.

### Example: Square with C<sub>4v</sub> Symmetry

For the square configuration space with 2 occupied and 2 unoccupied sites:

Without symmetry, there are $\binom{4}{2} = 6$ distinct configurations.

With C<sub>4v</sub> symmetry (8 operations: identity, 3 rotations, 4 reflections), these reduce to 2 unique configurations:

![Square unique configurations](./figures/square_unique_configurations.pdf)

1. **Adjacent sites occupied** - 4 equivalent configurations
2. **Diagonal sites occupied** - 2 equivalent configurations

Total: 2 unique configurations representing all 6 possible arrangements.

### Degeneracy

Each unique configuration has an associated **degeneracy** (or multiplicity) - the number of equivalent configurations it represents. This is stored in the `count` attribute:
```python
from bsym import ConfigurationSpace, PointGroup

# Create C4v symmetry group (8 operations)
c4v = PointGroup([e, c4, c4_inv, c2, sigma_x, sigma_y, sigma_ac, sigma_bd])

# Create configuration space
config_space = ConfigurationSpace(['a', 'b', 'c', 'd'], symmetry_group=c4v)

# Find unique configurations
unique_configs = config_space.unique_configurations({1: 2, 0: 2})

for config in unique_configs:
    print(f"{config.tolist()}: degeneracy = {config.count}")

# Output:
# [0, 0, 1, 1]: degeneracy = 4
# [0, 1, 0, 1]: degeneracy = 2
```

The degeneracies sum to the total number of configurations: 4 + 2 = 6.

## The Algorithm

### Overview

The bsym enumeration algorithm works by:

1. **Generate candidates**: Create all possible distinct permutations of the object labels
2. **Check for equivalence**: For each candidate, check if it's equivalent to any previously identified unique configuration
3. **Store unique configurations**: If not equivalent, add it to the set of unique configurations
4. **Track degeneracy**: Count how many configurations are equivalent to each unique configuration

### Key Algorithmic Insight

The algorithm maintains a "seen" set containing all configurations that have been identified (either as unique or as equivalent to a unique configuration).

For each candidate permutation:
- **If already in "seen"**: Skip it (it's equivalent to a previously found unique configuration)
- **If not in "seen"**: 
  1. This is a new unique configuration
  2. Generate all its symmetry equivalents by applying all symmetry operations
  3. Add all equivalents to the "seen" set
  4. Record the count of equivalents as the configuration's degeneracy

This ensures symmetry operations are only applied to configurations that represent genuinely new unique arrangements.

### Generating Candidates

For a given site distribution, the algorithm generates all unique permutations of object labels using Knuth's Algorithm L, which efficiently handles indistinguishable objects.

For example, with site distribution `{1: 2, 0: 2}`:
```
[0, 0, 1, 1]
[0, 1, 0, 1]
[0, 1, 1, 0]
[1, 0, 0, 1]
[1, 0, 1, 0]
[1, 1, 0, 0]
```

### Testing Equivalence with Hashing

Configurations are stored using their byte array representation, allowing fast hash-based lookup and comparison.

### Worked Example: Square with C<sub>2</sub> Rotation

Consider a simple square with only one non-identity symmetry operation: 180° rotation (C<sub>2</sub>), which maps positions 0↔2 and 1↔3.

We'll enumerate configurations with site distribution `{1: 2, 0: 2}`.

**Permutations to check**: `[0,0,1,1]`, `[0,1,0,1]`, `[0,1,1,0]`, `[1,0,0,1]`, `[1,0,1,0]`, `[1,1,0,0]`

---

**Step 1**: Process `[0, 0, 1, 1]`

1. Check if in "seen" set → No
2. This is a new unique configuration
3. Generate symmetry equivalents:
   - Identity: `[0, 0, 1, 1]`
   - C<sub>2</sub>: `[1, 1, 0, 0]` (positions 0↔2, 1↔3)
4. Add both to "seen" set: `{[0,0,1,1], [1,1,0,0]}`
5. Record unique configuration `[0, 0, 1, 1]` with degeneracy = 2

---

**Step 2**: Process `[0, 1, 0, 1]`

1. Check if in "seen" set → No
2. This is a new unique configuration
3. Generate symmetry equivalents:
   - Identity: `[0, 1, 0, 1]`
   - C<sub>2</sub>: `[0, 1, 0, 1]` (same! This configuration is symmetric)
4. Add to "seen" set: `{[0,0,1,1], [1,1,0,0], [0,1,0,1]}`
5. Record unique configuration `[0, 1, 0, 1]` with degeneracy = 1

---

**Step 3**: Process `[0, 1, 1, 0]`

1. Check if in "seen" set → No
2. This is a new unique configuration
3. Generate symmetry equivalents:
   - Identity: `[0, 1, 1, 0]`
   - C<sub>2</sub>: `[1, 0, 0, 1]`
4. Add both to "seen" set: `{[0,0,1,1], [1,1,0,0], [0,1,0,1], [0,1,1,0], [1,0,0,1]}`
5. Record unique configuration `[0, 1, 1, 0]` with degeneracy = 2

---

**Step 4**: Process `[1, 0, 0, 1]`

1. Check if in "seen" set → **Yes!**
2. Skip (already counted as equivalent to `[0, 1, 1, 0]`)

---

**Step 5**: Process `[1, 0, 1, 0]`

1. Check if in "seen" set → No
2. This is a new unique configuration
3. Generate symmetry equivalents:
   - Identity: `[1, 0, 1, 0]`
   - C<sub>2</sub>: `[1, 0, 1, 0]` (same! This configuration is symmetric)
4. Add to "seen" set
5. Record unique configuration `[1, 0, 1, 0]` with degeneracy = 1

---

**Step 6**: Process `[1, 1, 0, 0]`

1. Check if in "seen" set → **Yes!**
2. Skip (already counted as equivalent to `[0, 0, 1, 1]`)

---

**Result**: 4 unique configurations with degeneracies [2, 1, 2, 1]

Total configurations represented: 2 + 1 + 2 + 1 = 6 ✓

### Implementation in bsym

The same enumeration can be performed using bsym:
```python
from bsym import ConfigurationSpace, SymmetryGroup, SymmetryOperation

# Define C2 symmetry group
identity = SymmetryOperation.from_vector([1, 2, 3, 4], label='E')
c2 = SymmetryOperation.from_vector([3, 4, 1, 2], label='C2')
symmetry_group = SymmetryGroup([identity, c2])

# Create configuration space
config_space = ConfigurationSpace(
    objects=[1, 2, 3, 4],
    symmetry_group=symmetry_group
)

# Find unique configurations
unique_configs = config_space.unique_configurations({1: 2, 0: 2})

# Display results
for config in unique_configs:
    print(f"{config.tolist()}: degeneracy = {config.count}")

# Output:
# [0, 0, 1, 1]: degeneracy = 2
# [0, 1, 0, 1]: degeneracy = 1
# [0, 1, 1, 0]: degeneracy = 2
# [1, 0, 1, 0]: degeneracy = 1
```

The algorithm produces 4 unique configurations with the correct degeneracies, matching our manual enumeration.

## Site Distributions and Composition

### Specifying Object Distributions

The `unique_configurations()` method requires a **site distribution** dictionary that specifies how many objects of each type to place in the configuration space.

Format: `{label: count, ...}`

Examples:
- `{1: 2, 0: 2}` - 2 objects of type 1, 2 objects of type 0
- `{0: 1, 1: 2, 2: 1}` - ternary system with 1 of type 0, 2 of type 1, 1 of type 2
- `{1: 4}` - 4 objects of type 1 (remaining sites implicitly unoccupied)

### Label Conventions

Object labels are arbitrary integers. Common conventions:

- **Binary systems**: Use 0 and 1
  - 0 for vacant, 1 for occupied
  - 0 for species A, 1 for species B

- **Multi-species systems**: Use 0, 1, 2, 3, ...
  - Sequential integers for different species

The specific integers don't matter - what matters is which positions have the same vs different labels.

### Total Site Count

The sum of counts in the site distribution must equal the dimension of the configuration space:
```python
config_space = ConfigurationSpace(objects=[1, 2, 3, 4])  # 4-dimensional

# Valid:
config_space.unique_configurations({1: 2, 0: 2})  # 2 + 2 = 4 ✓

# Invalid:
config_space.unique_configurations({1: 3, 0: 2})  # 3 + 2 = 5 ✗
```

### Partial Occupancy

If you want some sites to remain unoccupied or unlabelled, include them explicitly:
```python
# 10-site system with 3 occupied sites and 7 vacant
config_space.unique_configurations({1: 3, 0: 7})
```

This is particularly relevant for crystallographic applications where not all sites are substituted.
