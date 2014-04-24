import numpy as np

def strip_first_col( fname, delimiter=None ):
    with open( fname, 'r' ) as fin:
        for line in fin:
            try:
               yield( line.split(delimiter, 1)[1] )
            except IndexError:
               continue

class SpaceGroup:
    def __init__( self, symmetry_operations = [] ):
        self.symmetry_operations = symmetry_operations

    @classmethod
    def read_from_file( cls, filename ):
        from symmetry import symmetry_operation as so
        data = np.loadtxt( filename, dtype=int )
        symmetry_operations = [ so.SymmetryOperation.from_vector( row.tolist() ) for row in data ]
        return( cls( symmetry_operations = symmetry_operations ) )

    @classmethod
    def read_from_file_with_labels( cls, filename ):
        from symmetry import symmetry_operation as so
        data = np.genfromtxt( filename, dtype=None )
        labels = [ row.tolist()[0].decode( 'utf-8' ) for row in data ]
        symmetry_operations = [ so.SymmetryOperation.from_vector( row.tolist()[1:] ) for row in data ]
        [ so.set_label( l ) for (l, so) in zip( labels, symmetry_operations ) ]
        return( cls( symmetry_operations = symmetry_operations ) )

    def save_symmetry_operation_vectors_to( self, filename ):
        operation_list = []
        for symmetry_operation in self.symmetry_operations:
            operation_list.append( symmetry_operation.as_vector() )
        np.savetxt( filename, np.array( operation_list ), fmt='%i' )

    def extend( self, symmetry_operations_list ):
        self.symmetry_operations.extend( symmetry_operations_list )
        return self

    def append( self, symmetry_operation ):
        self.symmetry_operations.append( symmetry_operation )
        return self

    def by_label( self, label ):
        return next((so for so in self.symmetry_operations if so.label == label), None)

    def labels( self ):
        return [ so.label for so in self.symmetry_operations ] 
