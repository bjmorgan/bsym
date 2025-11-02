# Configuration Spaces

## What is a Configuration Space?

A **configuration space** is an abstract vector space that represents possible arrangements of distinguishable objects. A configuration space consists of:

1. **A vector space** - n positions in a vector (indexed 0, 1, 2, ..., n-1)
2. **Object labels** - integers representing different types or species

This abstract framework can represent many different systems: atoms distributed across crystallographic sites, molecular conformations, states at discrete time steps, or any other system where discrete objects can be arranged in different ways.

### A Simple Example

Consider a 3-dimensional vector space. We can represent arrangements of two types of objects (labelled 0 and 1) as vectors like `[1, 1, 0]`.

If we interpret this vector as three sites in a triangle:

![Triangular configuration space](./figures/triangular_configuration_space.pdf)

then `[1, 1, 0]` represents:

![Triangular configuration example](./figures/triangular_configuration_example_1.pdf)

where positions 0 and 1 have object type 1 (shown in black), and position 2 has object type 0 (shown in white).

Different arrangements like `[0, 1, 1]` or `[1, 0, 1]` represent different **configurations** within the same **configuration space**.

### Configurations vs Configuration Spaces

- **Configuration Space**: The vector space itself (e.g., "3-dimensional space")
- **Configuration**: A specific assignment of labels (e.g., `[1, 1, 0]`)

## Mathematical Representation

### Configuration Vectors

Each configuration is represented as a vector of integers. For a configuration space with n positions, a configuration is an n-element vector:

$$\mathbf{v} = \begin{pmatrix}v_0\\v_1\\v_2\\\vdots\\v_{n-1}\end{pmatrix}$$

where each $v_i$ is a non-negative integer representing the type of object at position $i$.

### Object Labels

Object labels are arbitrary integers. Objects with the same label are considered indistinguishable. For example:

- Binary system: labels `0` and `1` (e.g., vacant/occupied, or species A/species B)
- Ternary system: labels `0`, `1`, and `2` (e.g., three different atomic species)
- Multi-species: any number of distinct integer labels

The specific integers used as labels are arbitrary - what matters is which positions have the same or different labels.

### Examples

For a 4-dimensional configuration space:

- `[0, 0, 1, 1]` - positions 0 and 1 have type 0, positions 2 and 3 have type 1
- `[0, 1, 0, 1]` - positions 0 and 2 have type 0, positions 1 and 3 have type 1
- `[2, 2, 1, 0]` - position 0 has type 2, position 1 has type 2, position 2 has type 1, position 3 has type 0
- `[1, 1, 1, 1]` - all positions have type 1

Each vector represents a distinct configuration. Whether two configurations like `[0, 0, 1, 1]` and `[0, 1, 0, 1]` are equivalent depends on the symmetry operations defined for the configuration space (discussed in the next section).

## The Configuration and ConfigurationSpace Classes

### The Configuration Class

In bsym, individual configurations are represented by `Configuration` objects. A `Configuration` stores:

- **A vector**: The integer array representing the configuration
- **Metadata**: Optional attributes like degeneracy counts

Creating a configuration:
```python
from bsym import Configuration

config = Configuration([1, 1, 0, 0])
```

### Numeric Representation

Each configuration has a numeric representation accessed via the `as_number` property. This provides a unique integer identifier for the configuration:
```python
config = Configuration([1, 2, 0])
print(config.as_number)  # Output: 120
```

This numeric representation is primarily used internally for efficient comparison and hashing of configurations during symmetry analysis.

### The ConfigurationSpace Class

A `ConfigurationSpace` object combines:

- **Objects**: A list defining the dimensionality of the space
- **Symmetry group**: Optional symmetry operations (defaults to identity only)

The objects list defines the vector space dimension:
```python
from bsym import ConfigurationSpace

# Create a 4-dimensional configuration space
config_space = ConfigurationSpace(objects=[1, 2, 3, 4])
```

The integers in the objects list serve as labels for the vector positions - they don't represent the configuration itself. They're often just sequential integers `[1, 2, 3, ..., n]`, but can be any distinct values.

### Configuration Space Without Symmetry

A `ConfigurationSpace` can be created without specifying symmetry operations. In this case, it contains only the identity operation, meaning no configurations are considered equivalent:
```python
config_space = ConfigurationSpace(objects=[1, 2, 3])
# Implicitly has only the identity symmetry operation
```

This is useful when you want to use the configuration space framework but don't need to identify symmetry-equivalent configurations.

## Why Use Abstract Representation?

### Separation of Concerns

The abstract vector representation separates the mathematical logic of symmetry analysis from the physical details of specific systems. This means:

- **Symmetry algorithms** work at the vector level, independent of coordinates or structures
- **Physical interpretation** is added as a separate layer when needed
- **The same code** handles crystals, molecules, or any other symmetric system

### Computational Efficiency

Working with integer vectors is computationally efficient:

- Integer comparisons are fast
- Vectors can be hashed and stored in sets/dictionaries
- No floating-point arithmetic or coordinate transformations needed during enumeration
- **Symmetry operations are simple permutations** of integer indices - just rearranging vector elements rather than matrix-vector multiplication with floating-point coordinates

### Generality

The abstract approach makes bsym applicable to any problem involving symmetric arrangements of discrete objects. You're not limited to crystallographic applications - the same framework handles:

- Disorder in crystal structures
- Molecular conformations
- Combinatorial problems with symmetry constraints
- Abstract group theory problems

### From Abstract to Physical

When working with real systems, the workflow is:

1. **Define the abstract configuration space** - vector dimension and symmetry operations
2. **Enumerate configurations** - find unique arrangements using vector-based algorithms
3. **Map to physical structures** - interpret abstract configurations as coordinates, structures, etc.

This separation allows the expensive symmetry analysis to happen at the abstract level, then efficiently generate corresponding physical structures only for the unique configurations.

## Connecting to Real Structures

### The CoordinateConfigSpace Class

For systems where vector positions correspond to physical coordinates, bsym provides `CoordinateConfigSpace`, which extends `ConfigurationSpace` with coordinate information:
```python
from bsym import CoordinateConfigSpace
import numpy as np

# Define coordinates for each position
coordinates = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])

# Create configuration space with coordinates
coord_space = CoordinateConfigSpace(coordinates, symmetry_group=my_symmetry_group)
```

The `CoordinateConfigSpace` maintains the abstract vector representation internally while also storing the associated coordinates. This allows symmetry analysis to happen at the abstract level, with results mapped back to coordinates when needed.

### The Pymatgen Interface

For crystallographic applications, bsym provides an interface to work with pymatgen `Structure` objects. This handles:

- Extracting symmetry operations from crystal structures
- Converting between abstract configurations and atomic structures
- Generating symmetry-inequivalent crystal structures from substitution patterns

The pymatgen interface is covered in detail in the [User Guide](../user_guide/index.rst). The key point is that it operates as a wrapper around the abstract `ConfigurationSpace` machinery - symmetry analysis happens at the vector level, then results are converted to `Structure` objects.
