import unittest
from bsym import bsym as sym

class IntegrationTest( unittest.TestCase ):

    def test_square_lattice_example( self ):
        spacegroup_filename = '../spacegroups/cubic_spacegroup_annotated'
        sg = sym.SpaceGroup.read_from_file( spacegroup_filename )
        site_dist = { 1 : 2, 
                      0 : 2 }
        unique_configurations = sym.process.unique_configurations_from_sites( site_dist, sg, verbose=False )
        self.assertEqual( len( unique_configurations ), 2 )
        c1, c2 = unique_configurations
        for nc1 in c1.numeric_equivalents( sg.symmetry_operations ):
            self.assertNotEqual( nc1, c2.tolist )
        for nc2 in c2.numeric_equivalents( sg.symmetry_operations ):
            self.assertNotEqual( nc1, c1.tolist )

if __name__ == '__main__':
    unittest.main()




