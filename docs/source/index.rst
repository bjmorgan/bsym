.. bsym documentation master file, created by
   sphinx-quickstart on Thu Jul  6 20:25:00 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

bsym - A basic symmetry module
==============================

|DOI| |Build Status| |Test Coverage| |Documentation Status|

``bsym`` is a basic Python symmetry module. It consists of core
classes that describe configuration vector spaces, their symmetry
operations, and specific configurations of objects withing these spaces.
The module also contains an interface for working with
`pymatgen <http://pymatgen.org>`__ ``Structure`` objects, to allow
simple generation of disordered symmetry-inequivalent structures from a
symmetric parent crystal structure.

API documentation is `here <modules.html>`__.

Examples are provided in a Jupyter notebook `here <http://nbviewer.jupyter.org/github/bjmorgan/bsym/blob/master/examples/bsym_examples.ipynb>`__

Source code is available as a git repository at https://github.com/bjmorgan/bsym.

Tests
-----

Automated testing of the latest commit happens
`here <https://travis-ci.org/bjmorgan/bsym>`__.

Manual tests can be run using

::

    python -m unittest discover

The code has been tested with Python versions 3.5 and above.

.. toctree::
   :caption: API documentation
   :maxdepth: 3

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |DOI| image:: https://zenodo.org/badge/19279643.svg
   :target: https://zenodo.org/badge/latestdoi/19279643
.. |Build Status| image:: https://travis-ci.org/bjmorgan/bsym.svg?branch=master
   :target: https://travis-ci.org/bjmorgan/bsym
.. |Test Coverage| image:: https://codeclimate.com/github/bjmorgan/bsym/badges/coverage.svg
   :target: https://codeclimate.com/github/bjmorgan/bsym/coverage
.. |Documentation Status| image:: https://readthedocs.org/projects/bsym/badge/?version=latest
   :target: http://bsym.readthedocs.io/en/latest/?badge=latest
