#! /usr/bin/env python3

from bsym import bsym as sym

sg_filename = "../../spacegroups/Li7Ti5O12_Li_surface"

sg = sym.SpaceGroup.read_from_file_with_labels( sg_filename )

c3_4 = sg.by_label( 'c3_4' )
sg.append( c3_4.invert().set_label( 'c3_4_i' ) )

s5_1 = sg.by_label( 's5_1' )
sg.append( s5_1.similarity_transform( c3_4 ).set_label( 's5_2' ) )
sg.append( s5_1.similarity_transform( c3_4.invert() ).set_label( 's5_2' ) )

# We have 2 Li ions to distribute amongst the [ 3 + 1 + 4 ] surface sites:
site_dist = { 1 : 2,
              0 : 6 }

configs = sym.process.unique_configurations_from_sites( site_dist, sg )
[ config.pprint() for config in configs ]

surface_sites_top = sym.SiteList.read_from_file( 'Li_top_sites' )
surface_sites_bottom = sym.SiteList.read_from_file( 'Li_bottom_sites' )

sym.process.coordinate_list_from_sitelists( configs = configs, 
                                            labels = [ 1 ], 
                                            sitelists = [ surface_sites_top, surface_sites_bottom ] ) 

