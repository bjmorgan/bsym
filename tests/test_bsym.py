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

    def test_bsym_imports_ColourOperation( self ):
        from bsym import ColourOperation

class TestOldBsymModule( unittest.TestCase ):

    def test_old_bsym_import_quits( self ):
        with self.assertRaises( SystemExit ):
            with warnings.catch_warnings():
                warnings.simplefilter( 'ignore' )
                from bsym import bsym

    def test_old_bsym_import_warns( self ):
        with self.assertRaises( Exception ) as w:
            with warnings.catch_warnings():
                warnings.simplefilter( 'error' )
                from bsym import bsym
                self.assertEqual( len(w), 1 )
                self.assertEqual( "You are trying to import bsym.bsym" in str( w[-1].message, True ) )

if __name__ == '__main__':
    unittest.main()
