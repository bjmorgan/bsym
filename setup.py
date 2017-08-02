import os
from bsym import __version__

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except ImportError:
    long_description = open('README.md').read()

config = {
    'description': 'A Basic Symmetry Module',
    'long_description': long_description,
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

setup(**config)
