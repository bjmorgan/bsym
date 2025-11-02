# Symmetry Operations

## What are Symmetry Operations?

A **symmetry operation** is a transformation of a configuration space that maps the vector positions onto each other. In the abstract vector representation, a symmetry operation is simply a permutation of the vector indices.

For a configuration space with n positions, a symmetry operation rearranges these positions while preserving the structure of the space. When applied to a configuration, a symmetry operation produces a new configuration by permuting the object labels according to the index permutation.

### Example: Triangular Configuration Space

Consider three positions arranged conceptually as an equilateral triangle:

![Triangular configuration space](./figures/triangular_configuration_space.pdf)

This is a 3-dimensional vector space with positions that we can label 0, 1, and 2 (corresponding to corners a, b, and c).

A specific configuration, such as `[1, 1, 0]`, assigns labels to each position:

![Triangular configuration example](./figures/triangular_configuration_example_1.pdf)

where positions 0 and 1 have object type 1 (shown in black), and position 2 has object type 0 (shown in white).

### Identity Operation

Every configuration space has at least one symmetry operation: the **identity operation** (often labelled E), which leaves all positions unchanged:

$$E = \begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{pmatrix}$$

This is the trivial permutation where every position maps to itself.

### Rotations and Reflections

For the triangular system, there are other symmetry operations including rotations and reflections:

![Triangular symmetry operations](./figures/triangular_example_symmetry_operations.pdf)

A reflection operation σ might swap positions, for example mapping position 1 to position 2 and position 2 to position 1, while leaving position 0 unchanged.

### Physical Interpretation

While symmetry operations are defined abstractly as permutations, they often correspond to physical transformations:

- **Rotations**: Rotating a structure around an axis
- **Reflections**: Mirroring across a plane
- **Inversions**: Mapping through a centre point
- **Translations**: Shifting by a lattice vector (for space groups)

These physical transformations, when applied to discrete sites, simply permute the site indices.

## Matrix and Vector Representation

### Permutation Matrices

A symmetry operation can be represented as a permutation matrix that transforms configuration vectors. For an n-dimensional space, this is an n × n matrix where each row and column contains exactly one 1, with all other elements 0.

For the identity operation on a 3-dimensional space:

$$E = \begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{pmatrix}$$

For a rotation operation that maps position 0→1, 1→2, 2→0:

$$C_3 = \begin{pmatrix}0 & 0 & 1\\1 & 0 & 0\\0 & 1 & 0\end{pmatrix}$$

Applying this matrix to configuration `[1, 1, 0]` gives:

$$\begin{pmatrix}0 & 0 & 1\\1 & 0 & 0\\0 & 1 & 0\end{pmatrix} \begin{pmatrix}1\\1\\0\end{pmatrix} = \begin{pmatrix}0\\1\\1\end{pmatrix}$$

### Vector Notation

A more compact representation uses a permutation vector showing where each position maps. Instead of the full matrix, we can write the C₃ rotation as:
```
[2, 3, 1]
```

meaning:
- Position 1 → Position 2
- Position 2 → Position 3
- Position 3 → Position 1

This vector notation is more compact and easier to read. In bsym, symmetry operations can be created from this vector form:
```python
from bsym import SymmetryOperation

# Create rotation operation using vector notation
c3 = SymmetryOperation.from_vector([2, 3, 1], label='C3')
```

Note that the vector notation uses 1-based indexing (positions numbered 1, 2, 3, ...), even though configurations use 0-based indexing internally (positions 0, 1, 2, ...).

### Applying Operations to Configurations

A symmetry operation transforms a configuration by permuting its elements:
```python
from bsym import Configuration, SymmetryOperation

config = Configuration([1, 1, 0])
c3 = SymmetryOperation.from_vector([2, 3, 1], label='C3')

# Apply the operation
result = c3.operate_on(config)
# Or equivalently:
result = c3 * config

print(result)  # Configuration([0, 1, 1])
```

![Rotation operation example](./figures/triangular_rotation_operation.pdf)

## Inverting Symmetry Operations

### Mathematical Inverse

For every symmetry operation A, there exists an inverse operation A⁻¹ such that:

$$A \cdot A^{-1} = E$$

where E is the identity operation. The inverse operation "undoes" the original operation.

### Example: C₃ and C₃⁻¹

For the triangular configuration space, the inverse of C₃ (clockwise rotation by 120°) is C₃⁻¹ (anticlockwise rotation by 120°):
```python
from bsym import SymmetryOperation

c3 = SymmetryOperation.from_vector([2, 3, 1], label='C3')
c3_inv = SymmetryOperation.from_vector([3, 1, 2], label='C3_inv')

# The product gives the identity
result = c3 * c3_inv
print(result.matrix)
# [[1. 0. 0.]
#  [0. 1. 0.]
#  [0. 0. 1.]]
```

![C3 inversion](./figures/triangular_c3_inversion.pdf)

### The `.invert()` Method

Symmetry operations can be inverted programmatically using the `.invert()` method:
```python
c3 = SymmetryOperation.from_vector([2, 3, 1], label='C3')

# Generate the inverse
c3_inv = c3.invert()
```

The resulting operation initially has no label. Labels can be set directly or using the `.set_label()` method:
```python
# Set label during inversion
c3_inv = c3.invert(label='C3_inv')

# Or set label afterwards
c3_inv = c3.invert().set_label('C3_inv')
```

## Symmetry Groups

### Definition

A **symmetry group** is a collection of symmetry operations. In bsym, a `SymmetryGroup` object contains a set of `SymmetryOperation` objects that describe the symmetries of a configuration space.
```python
from bsym import SymmetryGroup, SymmetryOperation

# Create individual operations
e = SymmetryOperation.from_vector([1, 2, 3], label='E')
c3 = SymmetryOperation.from_vector([2, 3, 1], label='C3')
c3_inv = SymmetryOperation.from_vector([3, 1, 2], label='C3_inv')

# Combine into a symmetry group
symmetry_group = SymmetryGroup([e, c3, c3_inv])
```

### Mathematical Groups

The term "symmetry group" is used by convention, but a `SymmetryGroup` in bsym is not required to be a complete mathematical group. It is simply a collection of symmetry operations that are relevant for a particular analysis.

A bsym `SymmetryGroup` may contain only a subset of operations - for example, you might include only the operations needed for a specific analysis, even if they don't form a complete group.

### The C<sub>3v</sub> Point Group Example

For the triangular configuration space, the complete C<sub>3v</sub> point group contains six operations: one identity, two rotations, and three reflections:

![C3v symmetry operations](./figures/triangular_c3v_symmetry_operations.pdf)
```python
from bsym import PointGroup, SymmetryOperation

# Identity
e = SymmetryOperation.from_vector([1, 2, 3], label='E')

# Rotations
c3 = SymmetryOperation.from_vector([2, 3, 1], label='C3')
c3_inv = SymmetryOperation.from_vector([3, 1, 2], label='C3_inv')

# Reflections
sigma_a = SymmetryOperation.from_vector([1, 3, 2], label='sigma_a')
sigma_b = SymmetryOperation.from_vector([3, 2, 1], label='sigma_b')
sigma_c = SymmetryOperation.from_vector([2, 1, 3], label='sigma_c')

# Create point group
c3v = PointGroup([e, c3, c3_inv, sigma_a, sigma_b, sigma_c])
```

### PointGroup and SpaceGroup

For convenience, bsym provides `PointGroup` and `SpaceGroup` classes that are functionally equivalent to `SymmetryGroup`:
```python
from bsym import PointGroup, SpaceGroup

# These are all equivalent:
group1 = SymmetryGroup([e, c3])
group2 = PointGroup([e, c3])
group3 = SpaceGroup([e, c3])
```

The different class names simply make code more readable by indicating what type of symmetry is being described.

## Creating Symmetry Groups from Files

### Reading Unlabelled Operations

Symmetry groups can be read from files containing permutation vectors. Each line represents one symmetry operation:
```
# example file: pair_symmetry.txt
1 2
2 1
```

This defines a symmetry group for two equivalent sites with the identity (1→1, 2→2) and a swap operation (1→2, 2→1).
```python
from bsym import SymmetryGroup

symmetry_group = SymmetryGroup.read_from_file('pair_symmetry.txt')
```

### Reading Labelled Operations

For better readability, operations can be labelled:
```
# example file: c3v_symmetry.txt
E       1 2 3
C3      2 3 1
C3_inv  3 1 2
sigma_a 1 3 2
sigma_b 3 2 1
sigma_c 2 1 3
```
```python
from bsym import SymmetryGroup

c3v = SymmetryGroup.read_from_file_with_labels('c3v_symmetry.txt')
```

The labels are used when displaying the symmetry group and can help identify specific operations during analysis.

### File Format

Both file formats use space-separated values:
- **Unlabelled**: Each line contains n integers representing the permutation
- **Labelled**: Each line starts with a label string, followed by n integers
- Comment lines starting with `#` are ignored
- Permutation vectors use 1-based indexing
