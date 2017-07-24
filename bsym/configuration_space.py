from bsym.permutations import flatten_list, unique_permutations
from bsym import Configuration, SymmetryGroup, SymmetryOperation
import numpy as np

class ConfigurationSpace:

    def __init__( self, objects, symmetry_group=None ):
        """
        Create a :any:`ConfigurationSpace` object.
  
        Args:
            objects (list): The set of objects that define the vector space of this configuration space.
            symmetry_group (:any:`SymmetryGroup`): The set of symmetry operations describing the symmetries of this configuration space.

        Returns:
            None
        """
        # Check that all properties have compatible dimensions
        dim = len( objects )
        self.objects = objects
        if symmetry_group:
            for so in symmetry_group.symmetry_operations:
                if so.matrix.shape[0] != dim:
                    raise ValueError
            self.symmetry_group = symmetry_group
        else:
            self.symmetry_group = SymmetryGroup( symmetry_operations=[ SymmetryOperation( np.identity( dim, dtype=int ), label='E' ) ] )

    def __repr__( self ):
        to_return = "ConfigurationSpace\n"
        to_return += self.objects.__repr__() + "\n"
        to_return += "\n".join( self.symmetry_group.__repr__().split("\n")[1:] )
        return to_return

    def unique_configurations( self, site_distribution, verbose=False ):
        """
        Find the symmetry inequivalent configurations for a given population of objects.

        Args:
            site_distribution (dict): A dictionary that defines the number of each object 
                                      to be arranged in this system.

                                      e.g. for a system with four sites, with two occupied (denoted `1`)
                                      and two unoccupied (denoted `0`)::

                                          { 1 : 2, 2 : 1 }
            verbose (opt:default=False): Print verbose output.

        Returns:
            unique_configurations (list): A list of :any:`Configuration` objects, for each symmetry 
                                          inequivalent configuration. 
        """
        if verbose:
            print( 'total number of sites: ' + str( sum( site_distribution.values() ) ) )
            print( 'using {:d} symmetry operations: '.format( len( self.symmetry_group.symmetry_operations ) ) )
        permutations = []
        working = True
        seen = set()
        unique_configurations = []
        s = flatten_list( [ [ key ] * site_distribution[ key ] for key in site_distribution ] )
        for new_permutation in unique_permutations( s ):
            if permutation_as_config_number( new_permutation) not in seen:
                config = Configuration.from_tuple( new_permutation )
                numeric_equivalents = set( config.numeric_equivalents( self.symmetry_group.symmetry_operations ) )
                config.count = len( numeric_equivalents )
                [ seen.add( i ) for i in numeric_equivalents ]
                unique_configurations.append( config )
                if verbose:
                    print( "found {:d}, screened {:d}".format( len( unique_configurations ), len( seen ) ) )
        if verbose:
            print( 'unique configurations: ' + str( len( unique_configurations ) ) )
        return( unique_configurations )

def permutation_as_config_number( p ):
    """
    A numeric representation of a numeric list.

    Example:
        >>> permutation_as_config_number( [ 1, 1, 0, 0, 1 ] )
        11001
    """
    tot = 0
    for num in p:
        tot *= 10
        tot += num
    return tot
