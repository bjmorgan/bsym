# Installation

## Requirements

bsym requires Python 3.10 or later.

## Standard Installation

Install from PyPI:
```bash
pip install bsym
```

## Installation from Source

Clone the repository and install:
```bash
git clone https://github.com/bjmorgan/bsym.git
cd bsym
pip install .
```

## Development Installation

For development work (running tests, building documentation):
```bash
git clone https://github.com/bjmorgan/bsym.git
cd bsym
pip install -e ".[dev]"
```

## Verifying Installation
```python
import bsym
print(bsym.__version__)
```
