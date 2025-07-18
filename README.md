# FuseFlow

A Python package for building FUSE hydrological models.

## Installation

### From PyPI (when published)
```bash
pip install fuseflow
```

After installation, the `fuseflow` command-line tool will be available in your environment.

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

### Python API

```python
from fuseflow import FUSEWorkflow

# Create a workflow instance
workflow = FUSEWorkflow(
    name="my_workflow",
    cat="/path/to/catchment.geojson",
    forcing_vars={"precip": "precipitation", "temp": "temperature"},
    forcing_files="/path/to/forcing/files/",
    forcing_units={"precip": "mm/day", "temp": "degC"},
    settings={"start_date": "2000-01-01", "end_date": "2010-12-31"}
)

# Run the workflow
result = workflow.run()
print(result)

# Save the results
workflow.save("/path/to/output/")
```

### Command Line Interface

The package also provides a CLI for running workflows from JSON configuration files:

```bash
# Run a workflow from a configuration file
fuseflow --json config.json --output-path ./results

# Run with verbose output
fuseflow --json config.json --output-path ./results --verbose

# Show help
fuseflow --help

# Show version
fuseflow --version
```

See [CLI_USAGE.md](CLI_USAGE.md) for detailed CLI documentation and examples.

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
