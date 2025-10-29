# bsym

[![PyPI version](https://badge.fury.io/py/bsym.svg)](https://badge.fury.io/py/bsym)
[![DOI](https://zenodo.org/badge/19279643.svg)](https://zenodo.org/badge/latestdoi/19279643)
[![status](http://joss.theoj.org/papers/6696543fc631bf66feb99a9cde808a39/status.svg)](http://joss.theoj.org/papers/6696543fc631bf66feb99a9cde808a39)
[![Coverage Status](https://coveralls.io/repos/github/bjmorgan/bsym/badge.svg?branch=master)](https://coveralls.io/github/bjmorgan/bsym?branch=master)
[![Documentation Status](https://readthedocs.org/projects/bsym/badge/?version=latest)](http://bsym.readthedocs.io/en/latest/?badge=latest)

`bsym` is a basic Python symmetry module. It consists of core classes that describe configuration vector spaces, their symmetry operations, and specific configurations of objects within these spaces. The module also contains an interface for working with [`pymatgen`](http://pymatgen.org) `Structure` objects, to allow simple generation of disordered symmetry-inequivalent structures from a symmetric parent crystal structure.

Usage examples are provided in the [documentation](http://bsym.readthedocs.io/en/latest).
API documentation is [here][API].
Source code is available at [https://github.com/bjmorgan/bsym][github].

## Requirements

`bsym` requires Python 3.10 or later.

## Installation

### Standard Installation

Install from PyPI:
```bash
pip install bsym
```

### Installation from Source

Download the latest release from [GitHub](https://github.com/bjmorgan/bsym/releases), and install:
```bash
cd bsym
pip install .
```

Or clone the latest development version:
```bash
git clone git@github.com:bjmorgan/bsym.git
cd bsym
pip install .
```

### Development Installation

To install with development dependencies (for running tests, type checking, building docs, etc.):
```bash
git clone git@github.com:bjmorgan/bsym.git
cd bsym
pip install -e ".[dev]"
```

This installs `bsym` in editable mode with additional tools for development.

## Tests

Tests use pytest. After installing with development dependencies, run:
```bash
pytest
```

For verbose output:
```bash
pytest -v
```

To run specific test files:
```bash
pytest tests/unit_tests/test_symmetry_group.py
```

## Example Usage

### Enumerating Symmetry-Inequivalent Structures

`bsym` can enumerate symmetry-inequivalent structures for disordered materials. Here's an example using a `pymatgen` `Structure` as input:
```python
from bsym.interface.pymatgen import unique_structure_substitutions
unique_structures = unique_structure_substitutions(
    parent_structure, 
    'X',  # Sites to substitute
    {'O': 8, 'F': 16}  # 8 oxygen, 16 fluorine
)
print(f"Found {len(unique_structures)} unique structures")
```

## Documentation

An overview of the capabilities of `bsym` along with example code is provided in the [documentation](http://bsym.readthedocs.io/en/latest/).

## Citing `bsym`

This code can be cited as:

Morgan, Benjamin J. (2017). *bsym - a Basic Symmetry Module*. The Journal of Open Source Software. http://doi.org/10.21105/joss.00370

### BibTeX
```bibtex
@article{Morgan_JOSS2017b,
  doi = {10.21105/joss.00370},
  url = {https://doi.org/10.21105/joss.00370},
  year  = {2017},
  month = {aug},
  publisher = {The Open Journal},
  volume = {2},
  number = {16},
  author = {Benjamin J. Morgan},
  title = {bsym: A basic symmetry module},
  journal = {The Journal of Open Source Software}
}
```

[github]: https://github.com/bjmorgan/bsym
[doi]: https://zenodo.org/badge/latestdoi/19279643
[API]: http://bsym.readthedocs.io/en/latest/modules.html
