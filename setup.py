import os
from bsym import __version__

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

config = {
    'description': 'A Basic Symmetry Module',
    'long_description': read('README.md'),
    'author': 'Benjamin J. Morgan',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/bsym',
    'download_url': "https://github.com/bjmorgan/bsym/archive/%s.tar.gz" % (__version__),
    'author_email': 'b.j.morgan@bath.ac.uk',
    'version': __version__,
    'install_requires': [ 'numpy', 
			  'pymatgen', 
                          'coverage',
                          'codeclimate-test-reporter' ],
    'license': 'MIT',
    'packages': [ 'bsym' ],
    'scripts': [],
    'name': 'bsym'
}
