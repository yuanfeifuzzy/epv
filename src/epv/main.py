import sys
import functools
from pathlib import Path
from datetime import datetime
from config import OPTIONS, COLUMNS
from collections.abc import Sequence

import polars as pl
from loguru import logger


def log_and_exit(message, code=1):
    logger.error(message)
    sys.exit(code)


def validate_data(table, **kwargs):
    options = {}
    table = Path(table)
    logger.debug(f'Loading data from {table} ...')
    df = pl.read_csv(table, separator='\t')
    logger.debug(f'Loaded {df.height:,} features from {table.name}')
    
    columns = df.columns
    missing_columns = set(COLUMNS) - set(columns)
    if missing_columns:
        logger.error(f'Table {table.name} has missing columns: {missing_columns}')
        sys.exit(1)
    
    count_columns, score_columns = kwargs.get('count_columns', []), kwargs.get('score_columns', [])
    if count_columns:
        if not set(count_columns).issubset(columns):
            ms = set(count_columns) - set(columns)
            log_and_exit(f'Specified count columns: {ms} not found in the table')
    else:
        count_columns = [column for column in columns if column.startswith('count_')]
        if not count_columns:
            log_and_exit('Failed to find count columns (columns start with count_)')
    count_columns = sorted(count_columns)
    
    if score_columns:
        if not set(score_columns).issubset(columns):
            ms = set(score_columns) - set(columns)
            log_and_exit(f'Specified zscore columns: {ms} not found in the table')
    else:
        score_columns = [column for column in columns if column.startswith('zscore_')]
        if not score_columns:
            log_and_exit('Failed to find zscore columns (columns start with zscore_')
    score_columns = sorted(score_columns)
    
    names = [column.replace('count_', '') for column in count_columns]
    n2 = [column.replace('zscore_', '') for column in score_columns]
    if names != n2:
        logger.error('Names of count and zscore columns mismatch:')
        log_and_exit(f'  {names} != {n2}')
        
    if len(names) < 2:
        log_and_exit(f'Insufficient count or zscore columns: {len(names)}, at least 2 columns needed')
    
    options['count_columns'] = count_columns
    options['score_columns'] = score_columns
    options['names'] = names
    
    
    x, y = kwargs.get('x', ''), kwargs.get('y', '')
    if x:
        if x not in names:
            log_and_exit(f'The specified x: {x} does not found in any count or zscore column')
    else:
        x = [name for name in names if name.upper() != 'NTC'][0]
        
    if y:
        if y not in names:
            log_and_exit(f'The specified y: {y} does not found in any count or zscore column')
    else:
        y = 'NTC' if 'NTC' in [name.upper() for name in names] else names[1]
    
    options['x'] = x
    options['y'] = y
    options['xc'], options['xs'] = f'count_{x}', f'zscore_{x}'
    options['yc'], options['ys'] = f'count_{y}', f'zscore_{y}'
    options['x_label'] = kwargs.get('x_label', f'zscore ({x})')
    options['y_label'] = kwargs.get('y_label', f'zscore ({y})')
    
    libraries = kwargs.get('libraries', [])
    libs = df['library'].unique().to_list()
    if libraries:
        ms = set(libraries) - set(libs)
        if ms:
            log_and_exit(f'Specified libraries: {ms} not found in table')
        options['libraries'] = list(libraries)
    else:
        options['libraries'] = libs
        
    colors = kwargs.get('colors', OPTIONS['colors']['default'])
    if not isinstance(colors, Sequence):
        log_and_exit(f'Invalid colors: {colors}. Colors must be a sequence of 3 hex color code strings')
    colors = list(colors)
    if len(colors) != 3:
        log_and_exit(f'Invalid colors: {colors}. Colors must be a sequence of 3 hex color code strings')
    options['colors'] = colors
    
    options['cards'] = kwargs.get('cards', OPTIONS['cards']['default'])
    options['date'] = kwargs.get('date', datetime.now().strftime('%m%d%Y'))
    
    options['outdir'] = Path(kwargs.get('outdir', table.parent.absolute()))
    if not options['outdir'].exists():
        options['outdir'].mkdir(exist_ok=True, parents=True)
    return df, options


def process_data(df, **kwargs):
    stat = {}
    xs, ys = kwargs['xs'], kwargs['ys']
    limits = df.group_by('library').agg(
            pl.col(xs).min().alias(f'x_min'),
            pl.col(xs).max().alias(f'x_max'),
            pl.col(ys).min().alias(f'y_min'),
            pl.col(ys).max().alias(f'y_max'),
    )
    stat['limits'] = {
        'min': min(limits['x_min'].min(), limits['y_min'].min()),
        'max': max(limits['x_max'].max(), limits['y_max'].max())
    }
    for row in limits.to_dicts():
        stat['limits'][row.pop('library')] = row


def overview_plot(df, **kwargs):
    image = kwargs['outdir'].joinpath('overview.png')
    if image.exists():
        logger.debug('Overview plot already exists')
    else:
        pass


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
    
    verbose = kwargs.get('verbose', OPTIONS['verbose']['default'])
    
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    formatter = '<level>[{time:YYYY-MM-DD HH:mm:ss}] {message}</level>'
    logger.add(sys.stdout, format=formatter, level=level)
    
    df, options = validate_data(table, **kwargs)
    process_data(df, **options)
    
    
if __name__ == '__main__':
    epv(Path('.').resolve().parent.parent / 'data/fake.data.tsv.gz')
    