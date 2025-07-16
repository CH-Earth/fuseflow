"""
This module provides functions to generate textual configuration files
using Jinja2 templating engine for FUSE model instantiations.
"""
# third-party libraries
from jinja2 import (
    Environment,
    PackageLoader,
)

from typing import (
    Dict,
    Sequence,
    Union
)

# define types
try:
    from os import PathLike
except ImportError:
    PathLike = str
else:
    PathLike = Union[str, PathLike]


# constants
TEMPLATE_SETTINGS = "fm_catch.txt.jinja"

# global variables and helper functions
def raise_helper(msg):
    """Jinja2 helper function to raise exceptions."""
    raise Exception(msg)
# Jinja2 environment setup
environment = Environment(
    loader=PackageLoader("fuseflow", "templates"),
    trim_blocks=True,
    lstrip_blocks=True,
    line_comment_prefix='##',
)
environment.globals['raise'] = raise_helper


def render_settings(
    paths: Dict[str, PathLike], # type: ignore
    dates: Dict[str, str],
    template_jinja_path: PathLike = TEMPLATE_SETTINGS
) -> str:
    """
    Render the settings file using the provided paths and dates.

    Parameters
    ----------
    paths : Dict[str, PathLike]
        A dictionary mapping path names to their file system locations.
    dates : Dict[str, str]
        A dictionary mapping date-related keys to their string values.
    template_jinja_path : PathLike, optional
        Path to the Jinja2 template file. Default is TEMPLATE_SETTINGS.

    Returns
    -------
    str
        The rendered settings file content as a string.
    """
    # create the template environment
    template = environment.get_template(template_jinja_path)

    # create content
    content = template.render(paths=paths, dates=dates)

    return content