# Introduction to bsym

## What is bsym?

bsym is a Python package for working with symmetry operations and enumerating symmetry-inequivalent configurations. It provides tools for:

- Defining abstract configuration spaces and their symmetry operations
- Finding all unique arrangements of objects that respect symmetry constraints
- Generating symmetry-inequivalent crystal structures with substitutional disorder

## What Problems Does bsym Solve?

### The Configuration Counting Problem

Consider a crystal with 16 anion sites where you want to substitute 8 oxygen atoms for 8 fluorine atoms. Without considering symmetry, there are $\binom{16}{8} = 12,870$ possible arrangements. However, if the crystal has symmetry operations (rotations, reflections, translations), many of these arrangements are equivalent - they're just the same structure viewed from different angles or shifted in space.

bsym identifies which arrangements are truly unique, dramatically reducing the number of structures you need to consider. For example, a high-symmetry crystal might have only a few dozen unique arrangements instead of thousands.

### Why This Matters

**For computational materials science:**
- Generate training sets for machine learning with only symmetry-inequivalent structures
- Calculate configuration-dependent properties efficiently
- Build phase diagrams by exploring composition space systematically
- Study defects and disorder without redundant calculations

**For mathematical and theoretical work:**
- Explore combinatorial problems with symmetry constraints
- Study group theory applications
- Develop enumeration algorithms

**For general research:**
- Avoid wasting computational resources on equivalent configurations
- Ensure systematic coverage of configuration space
- Calculate degeneracies for statistical mechanics

## Key Features

### Abstract Configuration Spaces

Work with symmetry at a mathematical level, independent of any physical system:
```python
from bsym import ConfigurationSpace, SymmetryGroup

config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
unique_configs = config_space.unique_configurations({1: 2, 0: 2})
```

### Crystallographic Interface

Integration with pymatgen for crystal structure generation:
```python
from bsym.interface.pymatgen import unique_structure_substitutions

unique_structures = unique_structure_substitutions(
    parent_structure, 'F', {'O': 8, 'F': 8}
)
```

### Efficient Algorithms

- Smart enumeration that only performs symmetry analysis when needed
- Species exchange optimization for composition enumeration
- Progress tracking for large systems

### Degeneracy Tracking

Each unique configuration includes its degeneracy - the number of symmetry-equivalent arrangements it represents. This is essential for statistical mechanics calculations.

## Who is bsym For?

**Computational materials scientists** studying:
- Solid solutions and substitutional disorder
- Defect configurations
- Surface adsorption patterns
- Magnetic ordering

**Researchers in related fields** working on:
- Combinatorial problems with symmetry
- Group theory applications
- Enumeration algorithms

**Students** learning about:
- Crystallographic symmetry
- Group theory
- Configuration spaces

## The Approach

bsym separates the mathematical logic of symmetry from the physical details of specific systems:

1. **Abstract representation**: Configurations are vectors of integers, symmetry operations are permutations
2. **Efficient computation**: Symmetry operations implemented as numpy array indexing, with hash-based configuration lookup
3. **Physical interpretation**: Map results back to structures, coordinates, etc. when needed

This separation makes the algorithms system-agnostic and computationally efficient.

## What's Next?

- **[Installation](installation.md)**: Get bsym installed
- **[Quickstart](quickstart.md)**: Try hands-on examples
- **[User Guide](../user_guide/index.rst)**: Detailed tutorials for specific tasks
- **[Theory](../theory/index.rst)**: Understand the concepts and algorithms

## Further Reading

For the theoretical background and detailed algorithm descriptions:
- [Configuration Spaces](../theory/configuration_spaces.md)
- [Symmetry Operations](../theory/symmetry_operations.md)
- [Unique Configuration Enumeration](../theory/unique_configurations.md)

For practical applications:
- [Basic Substitutions](../user_guide/basic_substitutions.ipynb)
- [Varying Composition](../user_guide/varying_composition.ipynb)
