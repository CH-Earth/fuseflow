# FuseFlow

A Python package for building FUSE hydrological models.

## Installation

### From PyPI (when published)
```bash
pip install fuseflow
```

### From Source
```bash
git clone https://github.com/kasra-keshavarz/fuseflow.git
cd fuseflow
pip install -e .
```

### Development Installation
```bash
git clone https://github.com/kasra-keshavarz/fuseflow.git
cd fuseflow
pip install -e ".[dev]"
```

## Quick Start

```python
from fuseflow import FuseFlow

# Create an instance
ff = FuseFlow()

# Use the package
result = ff.process()
print(result)
```

## Features

- Feature 1: Description
- Feature 2: Description
- Feature 3: Description

## Documentation

Full documentation is available at [link-to-docs].

## Development

### Setting up the development environment

1. Clone the repository:
   ```bash
   git clone https://github.com/kasra-keshavarz/fuseflow.git
   cd fuseflow
   ```

2. Install in development mode with dev dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

3. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

### Running Tests

```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=fuseflow --cov-report=html

# Run specific test file
pytest tests/test_core.py
```

### Code Quality

```bash
# Format code
black src tests

# Lint code
flake8 src tests

# Type checking
mypy src
```

### Building the Package

```bash
# Build source and wheel distributions
python -m build
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for your changes
5. Ensure all tests pass and code is properly formatted
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to the branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a list of changes.
