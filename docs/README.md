# Documentation

This directory contains the documentation for the FuseFlow package.

## Building Documentation

To build the documentation locally:

```bash
# Install documentation dependencies
pip install -e ".[docs]"

# Build the documentation
cd docs
make html

# View the documentation
open _build/html/index.html
```

## Documentation Structure

- `source/`: Source files for Sphinx documentation
- `_build/`: Generated documentation (not included in version control)
- `Makefile`: Make commands for building documentation

## Writing Documentation

We use Sphinx with reStructuredText format. Key files:

- `source/conf.py`: Sphinx configuration
- `source/index.rst`: Main documentation page
- `source/api.rst`: API reference
- `source/examples.rst`: Usage examples
