#! /usr/bin/env python3

from bsym import bsym as sym

sg_filename = "../../spacegroups/Li2TiO3_111"

sg = sym.SpaceGroup.read_from_file_with_labels( sg_filename )

T_1p = sg.by_label( 'T_1p' )
T_2p = sg.by_label( 'T_2p' )
T_p111 = sg.by_label( 'T_p111' )
T_m111 = sg.by_label( 'T_m111' )
i1 = sg.by_label( 'i' )

T_3p = ( T_1p * T_2p ).set_label( 'T_3p' )
sg.append( T_3p )

for so in [ T_1p, T_2p, T_3p ]:
    sg.append( so * T_p111 )
    sg.append( so * T_m111 )

gen = (so for so in sg.symmetry_operations.copy() if so.label not in  ['i', 'e'] )

[ sg.append( so * i1 ) for so in gen ]

# We have 4 Ti and 8 Li ions to distribute amongst the [ 12 ] sites:
site_dist = { 1 : 4,
              2 : 8 }

for config in sym.process.unique_configurations_from_sites( site_dist, sg ):
    config.pprint()
