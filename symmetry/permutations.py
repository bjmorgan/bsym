import itertools

def unique( sequence, idfun=None ): 
    # order preserving
    if idfun is None:
        def idfun( x ): return x
    seen = {}
    result = []
    for item in sequence:
        marker = idfun( item )
        if marker in seen: continue
        seen[ marker ] = 1
        result.append( item )
    return result

def unique_permutations( this_list ):
    return unique( itertools.permutations( this_list ) )

def flatten_list( this_list ):
    return [ item for sublist in this_list for item in sublist ]

def all_permutations( labels, number_of_sites=None ):
    # labels is a list of site occupations and their number in this system
    occupation_list = flatten_list( [ [ key ] * labels[key] for key in labels ] )
    if number_of_sites:
        assert( len( occupation_list ) == number_of_sites )
    return unique_permutations( occupation_list )
