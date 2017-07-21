import unittest
from unittest.mock import Mock

class TestBsymTopLevelClasses( unittest.TestCase ):

    def test_bsym_imports_SymmetryOperation( self ):
        from bsym import SymmetryOperation

    def test_bsym_imports_SymmetryGroup( self ):
        from bsym import SymmetryGroup

    def test_bsym_imports_SpaceGroup( self ):
        from bsym import SpaceGroup

    def test_bsym_imports_PointGroup( self ):
        from bsym import PointGroup

    def test_bsym_imports_ConfigurationSpace( self ):
        from bsym import ConfigurationSpace

if __name__ == '__main__':
    unittest.main()
