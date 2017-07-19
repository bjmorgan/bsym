import unittest
from unittest.mock import Mock, patch
from bsym import process
from bsym import SpaceGroup, SymmetryOperation, Configuration
from bsym.sitelist import SiteList
import numpy as np

class TestProcess( unittest.TestCase ):

    @patch( 'bsym.process.unique_permutations' )
    @patch( 'bsym.process.Configuration.from_tuple' )
    @patch( 'bsym.process.flatten_list' )
    def test_unique_configurations_from_sites( self, mock_flatten_list, mock_configuration_from_tuple, mock_unique_permutations ):
        mock_flatten_list.return_value = [ 1, 0 ]
        mock_unique_permutations.side_effect = [ [ [ 0, 1 ], [ 1, 0 ] ] ]
        s0 = Mock( spec=SymmetryOperation )
        s1 = Mock( spec=SymmetryOperation )
        site_distribution = { 1 : 1, 0 : 1 }
        spacegroup = SpaceGroup( symmetry_operations=[ s0, s1 ] )
        c1 = Mock( spec=Configuration )
        c1.numeric_equivalents.return_value = [ 1, 10 ]
        c1.as_number = 1
        c2 = Mock( spec=Configuration )
        c2.numeric_equivalents.return_value = [ 10, 1 ]
        c2.as_number = 10
        mock_configuration_from_tuple.side_effect = [ c1, c2 ]
        configs = process.unique_configurations_from_sites( site_distribution, spacegroup )
        self.assertEqual( len( configs ), 1 )
        self.assertEqual( configs[0], c1 )

    @patch( 'bsym.process.unique_permutations' )
    @patch( 'bsym.process.Configuration.from_tuple' )
    @patch( 'bsym.process.flatten_list' )
    def test_unique_configurations_from_sites_alt( self, mock_flatten_list, mock_configuration_from_tuple, mock_unique_permutations ):
        mock_flatten_list.return_value = [ 1, 1, 0 ]
        mock_unique_permutations.side_effect = [ [ [ 0, 1, 1 ], [ 1, 0, 1 ], [ 1, 1, 0 ] ] ]
        s0 = Mock( spec=SymmetryOperation )
        s1 = Mock( spec=SymmetryOperation )
        site_distribution = { 1 : 2, 0 : 1 }
        spacegroup = SpaceGroup( symmetry_operations=[ s0, s1 ] )
        c1 = Mock( spec=Configuration )
        c1.numeric_equivalents.return_value = [ 11, 11 ]
        c1.as_number = 11
        c2 = Mock( spec=Configuration )
        c2.numeric_equivalents.return_value = [ 101, 110 ]
        c2.as_number = 101
        c3 = Mock( spec=Configuration )
        c3.numeric_equivalents.return_value = [ 110, 101 ]
        c3.as_number = 110
        mock_configuration_from_tuple.side_effect = [ c1, c2, c3 ]
        configs = process.unique_configurations_from_sites( site_distribution, spacegroup )
        self.assertEqual( len( configs ), 2 )
        self.assertEqual( c1 in configs, True )
        self.assertEqual( c2 in configs, True )

    def test_cordinate_list_from_sitelists( self ):
        sitelist = SiteList( [ [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 1.0 ] ] )
        labels = [ 1, 2, 3 ]
        configs = [ Configuration( [[3], [2], [1] ] ) ]
        coords = process.list_of_coordinates_from_sitelists( configs, labels, [ sitelist ] )
        np.testing.assert_array_equal( coords[0], np.array( [ [ 0.0, 0.0, 1.0 ], [ 0.0, 1.0, 0.0 ], [ 1.0, 0.0, 0.0 ] ] ) )

if __name__ == '__main__':
    unittest.main() 
