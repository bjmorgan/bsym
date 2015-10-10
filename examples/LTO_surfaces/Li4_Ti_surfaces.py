#! /usr/bin/env python3

from bsym import bsym as sym

sg_filename = "../../spacegroups/Li4Ti5O12_Li_surface"

sg = sym.SpaceGroup.read_from_file_with_labels( sg_filename )

# We have [ 3 Li ] ions to distribute amongst the [ 6 + 2 + 2 + 6 ] surface sites:
site_dist = { 1 : 1,
              2 : 2,
              0 : 13 }

configs = sym.process.unique_configurations_from_sites( site_dist, sg, verbose=True )
[ config.pprint() for config in configs ]

sitelists = [ sym.SiteList.read_from_file( 'Li4_Ti_top_sites' ),
              sym.SiteList.read_from_file( 'Li4_Ti_bottom_sites' ) ]

sym.process.coordinate_list_from_sitelists( configs = configs,
                                            labels = [ 1, 2 ],
                                            sitelists = sitelists )

