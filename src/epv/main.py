import functools
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import polars as pl
from numpy.ma.core import maximum

from processor import process_data
from drawer import Card, position_cards, overview_plot, library_plot
from config import OPTIONS, logger, config_logger, log_and_exit, validate_data_and_options


def interface(func):
    """Decorator to inject OPTIONS into the docstring with uniform Sphinx formatting."""
    
    # Build the injected docstring in Sphinx :param: style
    kwargs_doc = "\nkwargs:"
    for name, config in OPTIONS.items():
        _type = config.get('type', str).__name__
        desc = config.get("help", "")
        default = config.get("default")
        
        # Format as :param <name>: <description> (Default: <val>)
        # Format as :type <name>: <type>
        kwargs_doc += f"\n    :param {name}: {desc}"
        if default is not None:
            kwargs_doc += f" (Default: {default})"
        kwargs_doc += f"\n    :type {name}: {_type}"

    # Append to the existing docstring
    if func.__doc__:
        func.__doc__ = func.__doc__.rstrip() + "\n" + kwargs_doc
    else:
        func.__doc__ = kwargs_doc
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        for key, value in kwargs.items():
            if key not in OPTIONS:
                raise TypeError(f"'{key}' is an invalid keyword argument for {func.__name__}")
            
            expected_type = OPTIONS[key].get("type", str)
            if expected_type and not isinstance(value, expected_type):
                type_name = getattr(expected_type, "__name__", str(expected_type))
                raise TypeError(f"Argument '{key}' must be {type_name}, got {type(value).__name__}")
        return func(*args, **kwargs)
    
    return wrapper


@interface
def epv(table: str | Path, **kwargs):
    """
    Enrichment Profile Visualizer.
    
    :param table: Path to a table containing enrichment data in TSV format.
    :type table: str | Path
    """
    
    config_logger(kwargs.get('verbose', OPTIONS['verbose']['default']))
    
    df, options = validate_data_and_options(table, **kwargs)
    unique = process_data(df, **options)
    # overview_plot(unique, **options)
    for library in options['libraries']:
        library_plot(df, unique, library, **options)
    
    
if __name__ == '__main__':
    # epv(Path('__file__').resolve().parent / 'data/fake.data.tsv.gz')
    epv(Path('__file__').resolve().parent / 'data/fake.data.tsv.gz', verbose=True, libraries=['qDOS36'])
    # position_cards(0, 10, [8, 9])
    