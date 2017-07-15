import unittest
from bsym import spacegroup
from bsym import symmetry_operation as so
from unittest.mock import Mock

class SpaceGroupTestCase( unittest.TestCase ):
    """Tests for SpaceGroup class"""

    def test_spacegroup_is_initialised( self ):
        s0, s1 = Mock( spec=so.SymmetryOperation ), Mock( spec=so.SymmetryOperation )
        sg = spacegroup.SpaceGroup( symmetry_operations=[ s0, s1 ] )
        self.assertEqual( sg.symmetry_operations[0], s0 )
        self.assertEqual( sg.symmetry_operations[1], s1 )

if __name__ == '__main__':
    unittest.main()
