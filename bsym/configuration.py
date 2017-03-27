import numpy as np

class Configuration( np.matrix ):

    def matches( self, test_configuration ):
        return ( ( self == test_configuration ).all() )

    def is_equivalent_to( self, test_configuration, symmetry_operations ):
        for symmetry_operation in symmetry_operations:
            if ( symmetry_operation.operate_on( self ).matches( test_configuration ) ):
                return True 
        else:
            return False

    def is_in_list( self, list ):
        for config in list:
            if ( self.matches( config ) ).all():
                return True
        else:
            return False

    def has_equivalent_in_list( self, list, symmetry_operations ):
        for config in list:
            print( config.as_number, self.as_number )
            if ( self.is_equivalent_to( config, symmetry_operations ) ):
                return True
        else:
            return False

    def set_lowest_numeric_representation( self, symmetry_operations ):
       self.lowest_numeric_representation = min( [ symmetry_operation.operate_on( self ).as_number for symmetry_operation in symmetry_operations ] )

    def numeric_equivalents( self, symmetry_operations ):
        return [ symmetry_operation.operate_on( self ).as_number for symmetry_operation in symmetry_operations ]

    @property
    def as_number( self ):
        return int( ''.join( str(e) for e in self.tolist() ) )

    @classmethod
    def from_tuple( cls, this_tuple ):
        return( cls( np.asarray( this_tuple ) ).T )

    def tolist( self ):
        return [ e[0] for e in super().tolist() ]

    def pprint( self ):
        print( ' '.join( [ str(e) for e in self.tolist() ] ) )

    def position( self, label ):
        return [ i for i,x in enumerate( self.tolist() ) if x == label ]
