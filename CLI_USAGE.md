# FuseFlow CLI Usage

The FuseFlow package includes a command-line interface (CLI) that allows you to run workflows from JSON configuration files.

## Installation

After installing the package, the CLI will be available as the `fuseflow` command.

## Basic Usage

```bash
fuseflow --json CONFIG_FILE --output-path OUTPUT_DIR [OPTIONS]
```

### Options

- `--json`, `-j`: Path to the JSON configuration file containing workflow parameters (required)
- `--output-path`, `-o`: Output directory path where the workflow results will be saved (required)
- `--verbose`, `-v`: Enable verbose output
- `--version`: Show the version and exit
- `--help`: Show help message and exit

## Configuration File Format

The configuration file should be a JSON file with the following structure:

```json
{
  "name": "example_workflow",
  "cat": "/path/to/catchment.geojson",
  "forcing_vars": {
    "precip": "precipitation",
    "temp": "temperature"
  },
  "forcing_files": "/path/to/forcing/files/",
  "forcing_units": {
    "precip": "mm/day",
    "temp": "degC"
  },
  "pet_method": "hamon",
  "model_number": [7, 37],
  "forcing_time_zone": "UTC",
  "model_time_zone": "America/Edmonton",
  "streamflow": "/path/to/streamflow.nc",
  "settings": {
    "start_date": "2000-01-01",
    "end_date": "2010-12-31",
    "start_evaluation": "2005-01-01",
    "end_evaluation": "2010-12-31"
  }
}
```

### Required Fields

- `cat`: Path to catchment GeoJSON file or GeoDataFrame
- `settings.start_date`: Start date for the simulation
- `settings.end_date`: End date for the simulation

### Optional Fields

- `name`: Name of the workflow (default: "fuse_models")
- `forcing_vars`: Dictionary mapping variable names to forcing file variables
- `forcing_files`: Path to forcing files directory or list of file paths
- `forcing_units`: Dictionary of units for forcing variables
- `pet_method`: Method for PET calculation (default: "hamon")
- `model_number`: Model number(s) to use (default: [7, 37])
- `forcing_time_zone`: Time zone for forcing data (default: "UTC")
- `model_time_zone`: Time zone for model (auto-detected if not provided)
- `streamflow`: Path to streamflow observation file
- `elev_bands`: Path to elevation bands file

## Examples

### Basic usage:
```bash
fuseflow --json config.json --output-path ./results
```

### With verbose output:
```bash
fuseflow --json config.json --output-path ./results --verbose
```

### Show help:
```bash
fuseflow --help
```

### Show version:
```bash
fuseflow --version
```

## Output Structure

The CLI will create the following directory structure in the output path:

```
output_path/
├── settings/          # Model configuration files
├── input/            # Input data files (NetCDF format)
├── output/           # Output directory for model results
├── file_manager_108.txt    # Model configuration for model 108
└── file_manager_2.txt   # Model configuration for model 2
```

## Error Handling

The CLI provides informative error messages for common issues:

- Missing or invalid configuration file
- Invalid JSON syntax
- Missing required configuration parameters
- Invalid parameter values
- File system errors

Use the `--verbose` flag to get detailed error information and stack traces for debugging.
