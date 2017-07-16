from bsym.permutations  import flatten_list, unique_permutations
from bsym.configuration import Configuration
import numpy as np

def unique_configurations_from_sites( site_distribution, spacegroup, verbose=False ):
    """
    Find the symmetry inequivalent configurations for a given population of objects and spacegroup.

    Args:
        site_distribution (dict): A dictionary that defines the number of each object 
                                  to be arranged in this system.

                                  e.g. for a system with four sites, with two occupied (denoted `1`)
                                  and two unoccupied (denoted `0`)::

                                      { 1 : 2, 2 : 1 }
        spacegroup (:any:`SpaceGroup`): The :any:`SpaceGroup` that desribes the symmetry 
                                        of this system.
        verbose (opt:default=False): Print verbose output.

    Returns:
        unique_configurations (list): A list of :any:`Configuration` objects, for each symmetry 
                                      inequivalent configuration. 
    """
    if verbose:
        print( 'total number of sites: ' + str( sum( site_distribution.values() ) ) )
        print( 'found {:d} inequivalent symmetry operations: '.format( len( spacegroup.symmetry_operations ) ) )
    permutations = []
    working = True
    seen = set()
    unique_configurations = []
    s = flatten_list( [ [ key ] * site_distribution[ key ] for key in site_distribution ] )
    for new_permutation in unique_permutations( s ):
        config = Configuration.from_tuple( new_permutation )
        if config.as_number not in seen:
            [ seen.add( i ) for i in config.numeric_equivalents( spacegroup.symmetry_operations ) ]
            unique_configurations.append( config )
            if verbose:
                print( "found {:d}, screened {:d}".format( len( unique_configurations ), len( seen ) ) )
    if verbose:
        print( 'unique configurations: ' + str( len( unique_configurations ) ) )
    return( unique_configurations ) 

def coordinate_list_from_sitelists( configs, labels, sitelists ):
    for idx, config in enumerate( configs, start=1 ):
        print( "\n# " + str( idx ) )
        for label in labels:
            for pos in config.position( label ):
                for sitelist in sitelists:
                    sitelist.print_site( pos )

def list_of_coordinates_from_sitelists( configs, labels, sitelists ):
    """
    Using
    1. A list of :any:`Configuration` objects,
    2. A list of labels,
    3. A list of :any:`SiteList` objects,
    for every :any:`Configuration`, select the positions by each label, in turn,
    and collect the corresponding coordinates from each :any:`SiteList`.
    
    Args:
        configs   (list): A list of :any:`Configuation` objects.
        labels    (list): A list of labels, e.g. [ 1, 2, 3 ].
        sitelists (list): A lsit of :any:`SiteList` objects.

    Returns:
        (np.array): A numpy array of coordinates.
    """ 
    all_coords = []
    for config in configs:
        coords = [] 
        for label in labels:
            for pos in config.position( label ):
                for sitelist in sitelists:
                    coords.append( sitelist[ pos ] )   
        all_coords.append( coords )
    return np.array( all_coords )
