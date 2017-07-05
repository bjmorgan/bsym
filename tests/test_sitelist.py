import unittest
from unittest.mock import Mock, mock_open, patch
import numpy as np
from bsym.sitelist import SiteList

import sys
from contextlib import contextmanager
from io import StringIO

@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err

class TestSiteList( unittest.TestCase ):

    def test_read_from_file( self ):
        coords = [ [ 1.0, 2.0, 3.0 ],
                   [ 4.0, 5.0, 6.0 ] ]
        with patch( 'numpy.loadtxt') as mock_loadtxt:
            mock_loadtxt.return_value = np.array( coords )
            sitelist = SiteList.read_from_file( 'filename' )
        self.assertEqual( type( sitelist ), SiteList )
        self.assertEqual( sitelist, coords )

    def test_print_site( self ):
        with captured_output() as ( out, err ):
            sitelist = SiteList( [ [ 1.0, 2.0, 3.0 ], [ 4.0, 5.0, 6.0 ] ] )
            sitelist.print_site(1)
        self.assertEqual( out.getvalue().strip(), '4.0 5.0 6.0' )

if __name__ == '__main__':
    unittest.main()
