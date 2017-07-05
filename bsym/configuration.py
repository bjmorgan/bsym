import numpy as np

class Configuration( np.matrix ):

    def matches( self, test_configuration ):
        """
        Test whether this configuration is equal to another configuration.

        Args:
            test_configuration (bsym.Configuration): The configuration to compare against.

        Returns:
            (bool): True | False.
        """
        return ( self == test_configuration ).all()

    def is_equivalent_to( self, test_configuration, symmetry_operations ):
        """
        Test whether this configuration is equivalent to another configuration
        under one or more of a set of symmetry operations.

        Args:
            test_configuration (bsym.Configuration): The configuration to compare against.
            symmetry_operations (list(bsym.SymmetryOperation): A list of SymmetryOperation objects.

        Returns:
            (bool): True | False
        """
        for symmetry_operation in symmetry_operations:
            if ( symmetry_operation.operate_on( self ).matches( test_configuration ) ):
                return True 
        else:
            return False

    def is_in_list( self, list ):
        """
        Test whether this configuration is in a list of configurations.

        Args:
            list (list(bsym.Comfiguration)): A list of Configuration instances.

        Returns:
            (bool): True | False
        """
        for config in list:
            if self.matches( config ):
                return True
        else:
            return False

    def has_equivalent_in_list( self, list, symmetry_operations ):
        """
        Test whether this configuration is equivalent by symmetry to one or more
        in a list of configurations.

        Args:
            list (list(bsym.Configuration)): A list of Configuration instances.
            symmetry_operations (list(bsym.SymmetryOperation)): A list of SymmetryOperation objects.

        Returns:
            (bool): True | False 
        """
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
