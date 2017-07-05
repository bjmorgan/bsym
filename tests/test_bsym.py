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

if __name__ == '__main__':
    unittest.main()
