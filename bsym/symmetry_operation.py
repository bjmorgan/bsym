from __future__ import annotations

from typing import overload
import numpy as np
from bsym.configuration import Configuration


def is_square(m: np.ndarray) -> bool:
    """
    Test whether a numpy matrix is square.

    Args:
        m: The matrix.

    Returns:
        True if matrix is square, False otherwise.
    """
    return bool(m.shape[0] == m.shape[1])


def is_permutation_matrix(m: np.ndarray) -> bool:
    """
    Test whether a numpy array is a `permutation matrix`_.

    .. _permutation_matrix: https://en.wikipedia.org/wiki/Permutation_matrix

    Args:
        m: The matrix.

    Returns:
        True if matrix is a permutation matrix, False otherwise.
    """
    m = np.asanyarray(m)
    return bool(
        m.ndim == 2 and m.shape[0] == m.shape[1] and
        (m.sum(axis=0) == 1).all() and
        (m.sum(axis=1) == 1).all() and
        ((m == 1) | (m == 0)).all()
    )


class SymmetryOperation:
    """
    `SymmetryOperation` class.
    """

    def __init__(
        self,
        matrix: np.ndarray | list[list[int]],
        label: str | None = None
    ) -> None:
        """
        Initialise a `SymmetryOperation` object

        Args:
            matrix: Square 2D vector as either a `numpy.ndarray` or `list`.
            label: Optional string label for this `SymmetryOperation` object.

        Raises:
            TypeError: if matrix is not `numpy.ndarray` or `list`.
            ValueError: if matrix is not square.
            ValueError: if matrix is not a `permutation matrix`_.

            .. _permutation_matrix: https://en.wikipedia.org/wiki/Permutation_matrix

        Notes:
            To construct a `SymmetryOperation` object from a vector of site mappings
            use the `SymmetryOperation.from_vector()` method.

        Returns:
            None
        """
        if not isinstance(matrix, (np.matrix, np.ndarray, list)):
            raise TypeError
        self.matrix = np.asarray(matrix)
        if not is_square(self.matrix):
            raise ValueError('Not a square matrix')
        if not is_permutation_matrix(self.matrix):
            raise ValueError('Not a permutation matrix')
        self.label = label
        self.index_mapping: np.ndarray = np.argmax(self.matrix, axis=1)

    @overload
    def __mul__(self, other: SymmetryOperation) -> SymmetryOperation: ...

    @overload
    def __mul__(self, other: Configuration) -> Configuration: ...

    def __mul__(self, other: SymmetryOperation | Configuration) -> SymmetryOperation | Configuration:
        """
        Multiply this `SymmetryOperation` matrix with another `SymmetryOperation`.

        Args:
            other: The other symmetry operation or configuration
                for the matrix multiplication self * other.

        Returns:
            A new `SymmetryOperation` instance with the resultant matrix,
            or a `Configuration` if `other` is a `Configuration`.
        """
        if isinstance(other, SymmetryOperation):
            return SymmetryOperation(self.matrix.dot(other.matrix))
        elif isinstance(other, Configuration):
            return self.operate_on(other)
        else:
            raise TypeError

    def invert(self, label: str | None = None) -> SymmetryOperation:
        """
        Invert this `SymmetryOperation` object.

        Args:
            label: Optional label for the inverted symmetry operation.

        Returns:
            A new `SymmetryOperation` object corresponding to the inverse matrix operation.
        """
        return SymmetryOperation(self.matrix.T.copy(), label=label)

    @classmethod
    def from_vector(
        cls,
        vector: list[int],
        count_from_zero: bool = False,
        label: str | None = None
    ) -> SymmetryOperation:
        """
        Initialise a SymmetryOperation object from a vector of site mappings.

        Args:
            vector: Vector of integers defining a symmetry operation mapping.
            count_from_zero: Set to True if the site index counts from zero.
            label: Optional string label for this `SymmetryOperation` object.

        Returns:
            A new SymmetryOperation object
        """
        if not count_from_zero:
            vector = [x - 1 for x in vector]
        dim = len(vector)
        matrix = np.zeros((dim, dim))
        for index, element in enumerate(vector):
            matrix[element, index] = 1
        new_symmetry_operation = cls(matrix, label=label)
        return new_symmetry_operation

    def similarity_transform(
        self,
        s: SymmetryOperation,
        label: str | None = None
    ) -> SymmetryOperation:
        """
        Generate the SymmetryOperation produced by a similarity transform S^{-1}.M.S

        Args:
            s: The symmetry operation or matrix S.
            label: The label to assign to the new SymmetryOperation.

        Returns:
            The SymmetryOperation produced by the similarity transform
        """
        s_new = s.invert() * (self * s)
        if label:
            s_new.set_label(label)
        return s_new

    def operate_on(self, configuration: Configuration) -> Configuration:
        """
        Return the Configuration generated by applying this symmetry operation

        Args:
            configuration: The configuration / occupation vector to operate on

        Returns:
            The new configuration obtained by operating on configuration with this 
            symmetry operation.
        """
        if not isinstance(configuration, Configuration):
            raise TypeError
        return Configuration(configuration.vector[self.index_mapping])

    def character(self) -> int:
        """
        Return the character of this symmetry operation (the trace of `self.matrix`).

        Args:
            None

        Returns:
            The trace of self.matrix
        """
        return int(np.trace(self.matrix))

    def as_vector(self, count_from_zero: bool = False) -> list[int]:
        """
        Return a vector representation of this symmetry operation

        Args:
            count_from_zero: Set to True if the vector representation counts from zero

        Returns:
            A vector representation of this symmetry operation (as a list)
        """
        offset = 0 if count_from_zero else 1
        return (np.argmax(self.matrix, axis=0) + offset).tolist()  # type: ignore[no-any-return]

    def set_label(self, label: str) -> SymmetryOperation:
        """
        Set the label for this symmetry operation.

        Args:
            label: Label to set for this symmetry operation

        Returns:
            self
        """
        self.label = label
        return self

    def pprint(self) -> None:
        """
        Pretty print for this symmetry operation

        Args:
            None

        Returns:
            None
        """
        label = self.label if self.label else '---'
        print(label + ' : ' + ' '.join([str(e) for e in self.as_vector()]))

    def __repr__(self) -> str:
        label = self.label if self.label else '---'
        return 'SymmetryOperation\nlabel(' + label + ")\n" + self.matrix.__repr__()