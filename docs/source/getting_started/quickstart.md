# Quickstart Guide

This guide provides quick examples to get you started with bsym, covering both abstract configuration spaces and practical crystallographic applications.

## Abstract Example: Symmetry-Inequivalent Arrangements

This example shows how to find unique arrangements of objects in a symmetric space without requiring any crystallographic knowledge.

### Problem: Four Sites in a Square

Consider four sites arranged in a square. How many unique ways can we place 2 occupied and 2 vacant sites, accounting for the square's rotational and reflection symmetry?

![Square configuration space](../theory/figures/square_configuration_space.pdf)

### Solution
```python
from bsym import ConfigurationSpace, SymmetryGroup, SymmetryOperation

# Define C4v symmetry operations (square symmetry)
e = SymmetryOperation.from_vector([1, 2, 3, 4], label='E')
c4 = SymmetryOperation.from_vector([2, 3, 4, 1], label='C4')
c4_inv = SymmetryOperation.from_vector([4, 1, 2, 3], label='C4i')
c2 = SymmetryOperation.from_vector([3, 4, 1, 2], label='C2')
sigma_x = SymmetryOperation.from_vector([4, 3, 2, 1], label='s_x')
sigma_y = SymmetryOperation.from_vector([2, 1, 4, 3], label='s_y')
sigma_ac = SymmetryOperation.from_vector([1, 4, 3, 2], label='s_ac')
sigma_bd = SymmetryOperation.from_vector([3, 2, 1, 4], label='s_bd')

# Create symmetry group
c4v = SymmetryGroup([e, c4, c4_inv, c2, sigma_x, sigma_y, sigma_ac, sigma_bd])

# Create configuration space
config_space = ConfigurationSpace(
    objects=['a', 'b', 'c', 'd'],
    symmetry_group=c4v
)

# Find unique configurations (2 occupied, 2 vacant)
unique_configs = config_space.unique_configurations({1: 2, 0: 2})

print(f"Found {len(unique_configs)} unique configurations")
for config in unique_configs:
    print(f"{config.tolist()}: degeneracy = {config.count}")
```

**Output:**
```
Found 2 unique configurations
[0, 0, 1, 1]: degeneracy = 4
[0, 1, 0, 1]: degeneracy = 2
```

Without symmetry, there would be 6 distinct arrangements. With C<sub>4v</sub> symmetry, these reduce to just 2 unique patterns:
- Adjacent sites occupied (4 equivalent arrangements)
- Diagonal sites occupied (2 equivalent arrangements)

### Next Steps

- Learn about [configuration spaces](../theory/configuration_spaces.md)
- Understand [symmetry operations](../theory/symmetry_operations.md)
- Read about the [enumeration algorithm](../theory/unique_configurations.md)

## Crystallographic Example: Disordered Structures

This example shows how to generate symmetry-inequivalent crystal structures with substitutional disorder.

### Problem: O/F Disorder in a Fluorite Structure

Generate all unique structures for a 2×2×1 CaF<sub>2</sub> supercell where we substitute 8 oxygen atoms for 8 of the 16 fluorine atoms.

### Solution
```python
from pymatgen.core import Structure
from bsym.interface.pymatgen import unique_structure_substitutions

# Load parent structure (CaF2 2x2x1 supercell)
parent_structure = Structure.from_file('CaF2_supercell.cif')
# Or create programmatically using pymatgen

# Generate all unique O/F arrangements
unique_structures = unique_structure_substitutions(
    structure=parent_structure,
    to_substitute='F',           # Substitute on F sites
    site_distribution={'O': 8, 'F': 8}  # 8 O, 8 F
)

print(f"Found {len(unique_structures)} symmetry-inequivalent structures")

# Check degeneracies
for i, structure in enumerate(unique_structures[:3]):
    n_equiv = structure.number_of_equivalent_configurations
    print(f"Structure {i}: represents {n_equiv} equivalent configurations")

# Export structures
for i, structure in enumerate(unique_structures):
    structure.to(filename=f'CaF2_O8F8_{i}.cif', fmt='cif')
```

**Output:**
```
Found 47 symmetry-inequivalent structures
Structure 0: represents 192 equivalent configurations
Structure 1: represents 192 equivalent configurations
Structure 2: represents 96 equivalent configurations
```

The `unique_structure_substitutions` function:
1. Automatically detects the space group symmetry
2. Identifies all symmetry-equivalent F sites
3. Enumerates only the symmetry-inequivalent O/F arrangements
4. Returns pymatgen `Structure` objects with degeneracy information

### Exploring Multiple Compositions

To generate structures across different O:F ratios:
```python
from bsym.interface.pymatgen import unique_structure_substitutions_by_composition

# Generate structures for all O:F compositions
all_structures = unique_structure_substitutions_by_composition(
    structure=parent_structure,
    to_substitute='F',
    species=['O', 'F']
)

# Results organized by composition
for composition, structures in all_structures.items():
    n_O, n_F = composition
    print(f"CaO{n_O}F{n_F}: {len(structures)} unique structures")
```

### Next Steps

- See [Basic Substitutions](../user_guide/basic_substitutions.ipynb) for more examples
- Learn about [varying composition](../user_guide/varying_composition.ipynb)
- Understand [degeneracy tracking](../user_guide/fixed_composition.ipynb)

## Key Concepts

### Configuration Space
An abstract vector space where you arrange different types of objects across discrete positions.

### Symmetry Operations
Transformations that map the configuration space onto itself (rotations, reflections, etc.).

### Symmetry-Inequivalent Configurations
The minimal set of configurations where no two can be transformed into each other by symmetry operations.

### Degeneracy
The number of symmetry-equivalent configurations represented by each unique configuration.

## Where to Go Next

**For abstract/mathematical applications:**
- [Theory: Configuration Spaces](../theory/configuration_spaces.md)
- [Theory: Symmetry Operations](../theory/symmetry_operations.md)

**For crystallographic applications:**
- [User Guide: Basic Substitutions](../user_guide/basic_substitutions.ipynb)
- [User Guide: Varying Composition](../user_guide/varying_composition.ipynb)

**For understanding the algorithms:**
- [Theory: Unique Configuration Enumeration](../theory/unique_configurations.md)
- [Theory: Composition Enumeration](../theory/composition_enumeration.md)
