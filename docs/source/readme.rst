|DOI| |Build Status| |Test Coverage| |Documentation Status|

``bsym`` is a basic Python symmetry module. It consists of some core
classes that describe configuration vector spaces, their symmetry
operations, and specific configurations of objects withing these spaces.
The module also contains an interface for working with
```pymatgen`` <http://pymatgen.org>`__ ``Structure`` objects, to allow
simple generation of disordered symmetry-inequivalent structures from a
symmetric parent crystal structure.

| API documentation is `here <http://bsym.readthedocs.io>`__.
| Examples are provided in a Jupyter notebook
  `here <http://nbviewer.jupyter.org/github/bjmorgan/bsym/blob/master/examples/bsym_examples.ipynb>`__
| Source code is available as a git repository at
  https://github.com/bjmorgan/bsym

Installation
------------

::

    pip install git+https://github.com/bjmorgan/bsym.git

Or download the latest release from
`GitHub <httpsL//github.com/bjmorgan/bsym/releases>`__, and install

::

    cd bsym
    python setup.py install

Or clone the latest development version

::

    git clone git@github.com:bjmorgan/bsym.git

and install the same way.

::

    cd bsym
    pythin setup.py install 

Documentation
-------------

An overview of the capabilities of ``bsym`` along with example code is
contained in a `Jupyter
notebook <http://jupyter-notebook.readthedocs.io/en/latest/#>`__ in the
repository ``examples`` directory
`examples/bsym\_examples.ipynb <http://nbviewer.jupyter.org/github/bjmorgan/bsym/blob/master/examples/bsym_examples.ipynb>`__.

API documentation is available `here <http://bsym.readthedocs.io>`__.

Tests
-----

Automated testing of the latest commit happens
`here <https://travis-ci.org/bjmorgan/bsym>`__.

Manual tests can be run using

::

    python -m unittest discover

The code has been tested with Python versions 3.5 and above.

Citing ``bsym``
---------------

This code can be cited as:

Morgan, Benjamin J. (2017, July 5). *bsym - a Basic Symmetry Module*.
Zenodo. http://doi.org/10.5281/zenodo.823127

BibTeX
~~~~~~

::

    @misc{morgan_benjamin_j_2017_823127,
      author       = {Morgan, Benjamin J.},
      title        = {bsym - a Basic Symmetry Module},
      month        = jul,
      year         = 2017,
      doi          = {10.5281/zenodo.823127},
      url          = {https://doi.org/10.5281/zenodo.823127}
    }

.. |DOI| image:: https://zenodo.org/badge/19279643.svg
   :target: https://zenodo.org/badge/latestdoi/19279643
.. |Build Status| image:: https://travis-ci.org/bjmorgan/bsym.svg?branch=master
   :target: https://travis-ci.org/bjmorgan/bsym
.. |Test Coverage| image:: https://codeclimate.com/github/bjmorgan/bsym/badges/coverage.svg
   :target: https://codeclimate.com/github/bjmorgan/bsym/coverage
.. |Documentation Status| image:: https://readthedocs.org/projects/bsym/badge/?version=latest
   :target: http://bsym.readthedocs.io/en/latest/?badge=latest
