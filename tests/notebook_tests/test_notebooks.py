import unittest
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import os
import glob

def get_cwd():
    path = os.path.abspath(__file__)
    return os.path.dirname(path)

def execute_notebook( filename ):
    with open( filename ) as f:
        nb = nbformat.read( f, as_version=4 )
        ep = ExecutePreprocessor( timeout=600, kernel_name='python3' )
        ep.preprocess( nb, { 'metadata': { 'path': os.path.dirname( filename ) } } )

def get_notebook_filenames( notebook_dir ):
    return glob.glob( os.path.join( get_cwd(), notebook_dir, '*.ipynb' ) )

class JupyterNotebookTestCase( unittest.TestCase ):

    def test_notebooks_execute( self ):
        notebook_dir = '../../examples'
        notebook_filenames = get_notebook_filenames( notebook_dir )
        for nf in notebook_filenames:
            execute_notebook( nf )

if __name__ == '__main__':
    unittest.main()
