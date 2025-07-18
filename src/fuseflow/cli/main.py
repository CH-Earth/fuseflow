"""Command line interface for FuseFlow."""

import sys
import json
import click
from pathlib import Path
from ..core import FUSEWorkflow


@click.command()
@click.option(
    '--json', 
    '-j', 
    'config_file',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Path to the JSON configuration file containing workflow parameters'
)
@click.option(
    '--output-path', 
    '-o', 
    type=click.Path(path_type=Path),
    required=True,
    help='Output directory path where the workflow results will be saved'
)
@click.option(
    '--verbose', 
    '-v', 
    is_flag=True, 
    help='Enable verbose output'
)
@click.version_option()
def cli(config_file: Path, output_path: Path, verbose: bool):
    """
    FuseFlow: A Python package for configuring FUSE hydrological models.
    
    Run a FUSE workflow from a JSON configuration file.
    """
    try:
        if verbose:
            click.echo(f"Loading configuration from: {config_file}")
        
        # Load the JSON configuration file
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        if verbose:
            click.echo(f"Configuration loaded successfully")
            click.echo(f"Creating FUSEWorkflow instance...")
        
        # Create the FUSEWorkflow instance from the configuration
        workflow = FUSEWorkflow.from_dict(config)
        
        if verbose:
            click.echo(f"FUSEWorkflow instance created: {workflow}")
            click.echo(f"Running workflow...")
        
        # Run the workflow
        result = workflow.run()
        
        if verbose:
            click.echo(f"Workflow completed: {result}")
            click.echo(f"Saving results to: {output_path}")
        
        # Save the results
        save_result = workflow.save(output_path)
        
        click.echo(f"✓ {save_result}")
        
    except FileNotFoundError as e:
        click.echo(f"Error: Configuration file not found: {e}", err=True)
        sys.exit(1)
    except json.JSONDecodeError as e:
        click.echo(f"Error: Invalid JSON in configuration file: {e}", err=True)
        sys.exit(1)
    except KeyError as e:
        click.echo(f"Error: Missing required configuration key: {e}", err=True)
        sys.exit(1)
    except ValueError as e:
        click.echo(f"Error: Invalid configuration value: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: An unexpected error occurred: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    cli()
