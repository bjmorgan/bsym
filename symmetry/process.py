from symmetry.permutations  import unique_permutations, all_permutations
from symmetry.configuration import Configuration

def unique_configurations_from_sites( site_distribution, spacegroup, verbose=False ):
    number_of_sites = sum( site_distribution.values() )
    if verbose:
        print( 'total number of sites: ' + str( number_of_sites ) ) 
    permutations = all_permutations( site_distribution, number_of_sites = number_of_sites )
    if verbose:
        print( 'total permutations: ' + str( len( permutations ) ) ) 
    all_configurations = [ Configuration.from_tuple( p ) for p in permutations ]
    unique_configurations = []
    for config in all_configurations:
        if not config.has_equivalent_in_list( unique_configurations, spacegroup.symmetry_operations ): 
            unique_configurations.append( config ) 
    if verbose:
        print( 'unique configurations: ' + str( len( unique_configurations ) ) ) 
    return( unique_configurations )

def coordinate_list_from_sitelists( configs, labels, sitelists ):
    for idx, config in enumerate( configs ):
        print( "\n# " + str( idx + 1 ) )
        for label in labels:
            for pos in config.position( label ):
                for sitelist in sitelists:
                    sitelist.print_site( pos )