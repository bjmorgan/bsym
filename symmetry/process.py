from symmetry.permutations  import unique_permutations, all_permutations
from symmetry.configuration import Configuration

def unique_configurations_from_sites( site_distribution, spacegroup ):
    number_of_sites = sum( site_distribution.values() )
    permutations = all_permutations( site_distribution, number_of_sites = number_of_sites )
    all_configurations = [ Configuration.from_tuple( p ) for p in permutations ]
    unique_configurations = []
    for config in all_configurations:
        if not config.has_equivalent_in_list( unique_configurations, spacegroup.symmetry_operations ): 
            unique_configurations.append( config ) 
    return( unique_configurations )