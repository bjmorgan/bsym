import numpy as np
from bsym import SymmetryOperation 

class SymmetryGroup:
    """
    :any:`SymmetryGroup` class.

    A :any:`SymmetryGroup` object contains a set of :any:`SymmetryOperation` objects.

    e.g.::

        SymmetryGroup( symmetry_operations=[ s1, s2, s3 ] )

    where `s1`, `s2`, and `s3` are :any:`SymmetryOperation` objects.

    :any:`SymmetryGroup` objects can also be created from files using the class methods::

        SymmetryGroup.read_from_file( filename )

    and::

        SymmetryGroup.read_from_file_with_labels( filename )
    """

    class_str = 'SymmetryGroup'

    def __init__( self, symmetry_operations=[] ):
        """
        Create a :any:`SymmetryGroup` object.

        Args:
            symmetry_operations (list): A list of :any:`SymmetryOperation` objects.

        Returns:
            None
        """
        self.symmetry_operations = symmetry_operations

    @classmethod
    def read_from_file( cls, filename ):
        """
        Create a :any:`SymmetryGroup` object from a file.
       
        The file format should be a series of numerical mappings representing each symmetry operation.

        e.g. for a pair of equivalent sites::

            # example input file to define the spacegroup for a pair of equivalent sites
            1 2
            2 1

        Args:
            filename (str): Name of the file to be read in.

        Returns:
            spacegroup (SymmetryGroup)	
        """
        data = np.loadtxt( filename, dtype=int )
        symmetry_operations = [ SymmetryOperation.from_vector( row.tolist() ) for row in data ]
        return( cls( symmetry_operations = symmetry_operations ) )

    @classmethod
    def read_from_file_with_labels( cls, filename ):
        """
        Create a :any:`SymmetryGroup` object from a file, with labelled symmetry operations.

        The file format should be a series of numerical mappings representing each symmetry operation, prepended with a string that will be used as a label.

        e.g. for a pair of equivalent sites::

            # example input file to define the spacegroup for a pair of equivalent sites
            E  1 2
            C2 2 1

        Args:
            filename (str): Name of the file to be read in.

        Returns:
            spacegroup (SymmetryGroup)
        """
        data = np.genfromtxt( filename, dtype=str )
        labels = [ row[0] for row in data ]
        vectors = [ [ float(s) for s in row[1:] ] for row in data ]
        symmetry_operations = [ SymmetryOperation.from_vector( v ) for v in vectors ]
        [ so.set_label( l ) for (l, so) in zip( labels, symmetry_operations ) ]
        return( cls( symmetry_operations=symmetry_operations ) )

    def save_symmetry_operation_vectors_to( self, filename ):
        """
        Save the set of vectors describing each symmetry operation in this :any:`SymmetryGroup` to a file.

        Args:
            filename (str): Name of the file to save to.

        Returns:
            None
        """
        operation_list = []
        for symmetry_operation in self.symmetry_operations:
            operation_list.append( symmetry_operation.as_vector() )
        np.savetxt( filename, np.array( operation_list ), fmt='%i' )

    def extend( self, symmetry_operations_list ):
        """
        Extend the list of symmetry operations in this :any:`SymmetryGroup`.

        Args:
            symmetry_operations_list (list): A list of :any:`SymmetryOperation` objects.

        Returns:
            self (:any:`SymmetryGroup`)
        """
        self.symmetry_operations.extend( symmetry_operations_list )
        return self

    def append( self, symmetry_operation ):
        """
        Append a :any:`SymmetryOperation` to this :any:`SymmetryGroup`.

        Args:
            symmetry_operation (:any:`SymmetryOperation`): The :any:`SymmetryOperation` to add.

        Returns:
            self (:any:`SymmetryGroup`)
        """
        self.symmetry_operations.append( symmetry_operation )
        return self

    def by_label( self, label ):
        """
        Returns the :any:`SymmetryOperation` with a matching label.

        Args:
            label (str): The label identifying the chosen symmetry operation.

        Returns:
            (:any:`SymmetryOperation`): The symmetry operation that matches this label.
        """ 
        return next((so for so in self.symmetry_operations if so.label == label), None)

    @property
    def labels( self ):
        """
        A list of labels for each :any:`SymmetryOperation` in this spacegroup.

        Args:
            None

        Returns:
            (list): A list of label strings.
        """
        return [ so.label for so in self.symmetry_operations ] 

    def __repr__( self ):
        to_return = '{}\n'.format( self.__class__.class_str )
        for so in self.symmetry_operations:
            to_return += "{}\t{}\n".format( so.label, so.as_vector() )
        return to_return

    @property
    def size( self ):
        return len( self.symmetry_operations )

