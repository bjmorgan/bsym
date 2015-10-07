#! /usr/bin/env python3

from bsym import bsym as sym

def complete_spacegroup():
    sg_filename = "../spacegroups/Li7Ti5O12_Li_surface"
    sg = sym.SpaceGroup.read_from_file_with_labels( sg_filename )
    # generate symmetry operations that aren't already in the sg file
    c3_4 = sg.by_label( 'c3_4' )
    sg.append( c3_4.invert().set_label( 'c3_4_i' ) )
    s5_1 = sg.by_label( 's5_1' )
    sg.append( s5_1.similarity_transform( c3_4 ).set_label( 's5_2' ) )
    sg.append( s5_1.similarity_transform( c3_4.invert() ).set_label( 's5_2' ) )
    return sg

# We have [ 1 Li + 1 Ti ] ions to distribute amongst the [ 3 + 1 + 4 ] surface sites:
site_dist = { 1 : 1,
              2 : 1,
              0 : 6 }

configs = sym.process.unique_configurations_from_sites( site_dist, complete_spacegroup() )
# [ config.pprint() for config in configs ]

sitelists = [ sym.SiteList.read_from_file( 'LiTi_top_sites' ),
              sym.SiteList.read_from_file( 'LiTi_bottom_sites' ) ]

sym.process.coordinate_list_from_sitelists( configs = configs, 
                                            labels = [ 1, 2 ], 
                                            sitelists = sitelists )         
