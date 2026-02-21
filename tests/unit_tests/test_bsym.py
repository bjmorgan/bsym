import unittest
from unittest.mock import Mock
import warnings 

class TestBsymTopLevelClasses( unittest.TestCase ):

    def test_bsym_imports_SymmetryOperation( self ):
        from bsym import SymmetryOperation

    def test_bsym_imports_SymmetryGroup( self ):
        from bsym import SymmetryGroup

    def test_bsym_imports_SpaceGroup( self ):
        from bsym import SpaceGroup

    def test_bsym_imports_PointGroup( self ):
        from bsym import PointGroup

    def test_bsym_imports_Configuration( self ):
        from bsym import Configuration

    def test_bsym_imports_ConfigurationSpace( self ):
        from bsym import ConfigurationSpace

    def test_bsym_imports_CoordinateConfigSpace( self ):
        from bsym import CoordinateConfigSpace

class TestOldBsymModule( unittest.TestCase ):

    def test_old_bsym_import_raises_import_error( self ):
        import sys
        # Remove cached module if present so the import re-executes
        sys.modules.pop('bsym.bsym', None)
        with self.assertRaises( ImportError ) as cm:
            from bsym import bsym
        self.assertIn("bsym.bsym", str(cm.exception))

if __name__ == '__main__':
    unittest.main()
