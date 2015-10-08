from bsym.permutations  import unique_permutations, all_permutations
from bsym.configuration import Configuration

def unique_configurations_from_sites( site_distribution, spacegroup, verbose=False ):
    number_of_sites = sum( site_distribution.values() )
    if verbose:
        print( 'total number of sites: ' + str( number_of_sites ) ) 
    permutations = all_permutations( site_distribution, number_of_sites = number_of_sites )
    if verbose:
        print( 'total permutations: ' + str( len( permutations ) ) ) 
    all_configurations = [ Configuration.from_tuple( p ) for p in permutations ]
    for i, config in enumerate( all_configurations ):
        config.set_lowest_numeric_representation( spacegroup.symmetry_operations ) 
    seen = set()
    unique_configuration_counts = {}
    unique_configurations = []
    for config in all_configurations:
        if config.lowest_numeric_representation not in seen:
            seen.add( config.lowest_numeric_representation )
            unique_configurations.append( config )
            unique_configuration_counts[ config.lowest_numeric_representation ] = 1
        else:
            unique_configuration_counts[ config.lowest_numeric_representation ] += 1
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
