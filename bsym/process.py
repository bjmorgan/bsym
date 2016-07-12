from bsym.permutations  import unique_permutations, next_permutationS, flatten_list
from bsym.configuration import Configuration
import numpy as np

def unique_configurations_from_sites( site_distribution, spacegroup, verbose=False ):
    number_of_sites = sum( site_distribution.values() )
    if verbose:
        print( 'total number of sites: ' + str( number_of_sites ) )
        print( 'found {:d} inequivalent symmetry operations: '.format( len( spacegroup.symmetry_operations ) ) )
    permutations = []
    working = True
    seen = set()
    unique_configurations = []
    count = 0
    s = flatten_list( [ [ key ] * site_distribution[ key ] for key in site_distribution ] )
    while working:
        working, new_permutation = next_permutationS( s )
        count += 1
        config = Configuration.from_tuple( new_permutation )
        config_id = int( ''.join( map( str, new_permutation ) ) )
        if config_id not in seen:
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
    all_coords = []
    for config in configs:
        coords = [] 
        for label in labels:
            for pos in config.position( label ):
                for sitelist in sitelists:
                    coords.append( sitelist[ pos ] )   
        all_coords.append( coords )
    return np.array( all_coords )
