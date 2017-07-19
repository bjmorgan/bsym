import unittest
from unittest.mock import Mock

class TestBsymImports( unittest.TestCase ):

    def test_bsym_imports_SpaceGroup( self ):
        from bsym.bsym import SpaceGroup

    def test_bsym_imports_Configuration( self ):
        from bsym.bsym import Configuration

    def test_bsym_imports_SiteList( self ):
        from bsym.bsym import SiteList

    def test_bsym_imports_process( self ):
        from bsym.bsym import process

    def test_bsym_imports_SymmetryOperation( self ):
        from bsym.bsym import SymmetryOperation

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
