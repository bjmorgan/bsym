import numpy as np
from bsym.configuration import Configuration

def is_square( m ):
    """
    Test whether a numpy matrix is square.

    Args:
        m (np.matrix): The matrix.

    Returns:
        (bool): True | False.
    """
    return m.shape[0] == m.shape[1]

def is_permutation_matrix( m ):
    """
    Test whether a numpy array is a `permutation matrix`_.

    .. _permutation_matrix: https://en.wikipedia.org/wiki/Permutation_matrix
    
    Args:
        m (mp.matrix): The matrix.

    Returns:
        (bool): True | False.
    """
    m = np.asanyarray(m)
    return (m.ndim == 2 and m.shape[0] == m.shape[1] and
            (m.sum(axis=0) == 1).all() and 
            (m.sum(axis=1) == 1).all() and
            ((m == 1) | (m == 0)).all())
    
class SymmetryOperation:
    """
    `SymmetryOperation` class.
    """

    def __init__( self, matrix, label=None ):
        """
        Initialise a `SymmetryOperation` object

        Args:
            matrix (numpy.matrix|numpy.ndarray|list): square 2D vector as either a
            `numpy.matrix`, `numpy.ndarray`, or `list`.
            for this symmetry operation.
            label (default=None) (str): optional string label for this `SymmetryOperation` object.
        Raises:
            TypeError: if matrix is not `numpy.matrix`, `numpy.ndarray`, or `list`.
            ValueError: if matrix is not square.
            ValueError: if matrix is not a `permutation matrix`_.

            .. _permutation_matrix: https://en.wikipedia.org/wiki/Permutation_matrix

        Notes:
            To construct a `SymmetryOperation` object from a vector of site mappings
            use the `SymmetryOperation.from_vector()` method.

        Returns:
            None
        """
        if isinstance( matrix, np.matrix ):
            self.matrix = np.array( matrix )
        elif isinstance( matrix, np.ndarray ):
            self.matrix = np.array( matrix )
        elif isinstance( matrix, list):
            self.matrix = np.array( matrix )
        else:
            raise TypeError
        if not is_square( self.matrix ):
            raise ValueError('Not a square matrix')
        if not is_permutation_matrix( self.matrix ):
            raise ValueError('Not a permutation matrix')
        self.label = label
        self.index_mapping = np.array( [ np.array(row).tolist().index(1) for row in matrix ] )

    def __mul__( self, other ):
        """
        Multiply this `SymmetryOperation` matrix with another `SymmetryOperation`.

        Args:
            other (SymmetryOperation, Configuration): the other symmetry operation or configuration or matrix
            for the matrix multiplication self * other.

        Returns:
            (SymmetryOperation): a new `SymmetryOperation` instance with the resultant matrix.
            (Configuration): if `other` is a `Configuration`.
        """
        if isinstance( other, SymmetryOperation ):
            return SymmetryOperation( self.matrix.dot( other.matrix ) )
        elif isinstance( other, Configuration ):
            return self.operate_on( other )
        else:
            raise TypeError

    def invert( self, label=None ):
        """
        Invert this `SymmetryOperation` object.

        Args:
            None
 
        Returns:
            A new `SymmetryOperation` object corresponding to the inverse matrix operation.
        """
        return SymmetryOperation( np.linalg.inv( self.matrix ).astype( int ), label=label )

    @classmethod
    def from_vector( cls, vector, count_from_zero=False, label=None ):
        """
        Initialise a SymmetryOperation object from a vector of site mappings.

        Args:
            vector (list): vector of integers defining a symmetry operation mapping.
            count_from_zero (default = False) (bool): set to True if the site index counts from zero.
            label (default=None) (str): optional string label for this `SymmetryOperation` object.
   
        Returns:
            a new SymmetryOperation object
        """
        if not count_from_zero:
            vector = [ x - 1 for x in vector ]
        dim = len( vector )
        matrix = np.zeros( ( dim, dim ) )
        for index, element in enumerate( vector ):
            matrix[ element, index ] = 1
        new_symmetry_operation = cls( matrix, label=label )
        return new_symmetry_operation

    def similarity_transform( self, s, label=None ):
        """
        Generate the SymmetryOperation produced by a similarity transform S^{-1}.M.S

        Args:
            s: the symmetry operation or matrix S.
            label (:obj:`str`, optional): the label to assign to the new SymmetryOperation. Defaults to None.

        Returns:
            the SymmetryOperation produced by the similarity transform
        """
        s_new = s.invert() * ( self * s )
        if label:
            s_new.set_label( label )
        return s_new

    def operate_on( self, configuration ):
        """
        Return the Configuration generated by appliying this symmetry operation

        Args:
            configuration (Configuration): the configuration / occupation vector to operate on

        Returns:
            (Configuration): the new configuration obtained by operating on configuration with this symmetry operation. 
        """
        if not isinstance( configuration, Configuration ):
            raise TypeError
        return Configuration( configuration.vector[ self.index_mapping ] )
        #return Configuration( self.matrix.dot( configuration.vector ) )

    def character( self ):
        """
        Return the character of this symmetry operation (the trace of `self.matrix`).

        Args:
            none

        Returns:
            np.trace( self.matrix )
        """
        return np.trace( self.matrix )

    def as_vector( self, count_from_zero=False ):
        """
        Return a vector representation of this symmetry operation

        Args:
            count_from_zero (default = False) (bool): set to True if the vector representation counts from zero
      
        Returns:
            a vector representation of this symmetry operation (as a list)
        """
        offset = 0 if count_from_zero else 1
        return [ row.tolist().index( 1 ) + offset for row in self.matrix.T ]

    def set_label( self, label ):
        """
        Set the label for this symmetry operation.
  
        Args:
            label: label to set for this symmetry operation
        Returns:
            self 
        """
        self.label = label
        return self

    def pprint( self ):
        """
        Pretty print for this symmetry operation

        Args:
            None
        Returns:
            None
        """
        label = self.label if self.label else '---'
        print( label + ' : ' + ' '.join( [ str(e) for e in self.as_vector() ] ) )
        
    def __repr__( self ):
        label = self.label if self.label else '---'
        return 'SymmetryOperation\nlabel(' + label + ")\n" + self.matrix.__repr__()
