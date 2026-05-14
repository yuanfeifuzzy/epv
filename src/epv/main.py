import sys
import functools
from pathlib import Path
from datetime import datetime
from config import OPTIONS, COLUMNS
from collections.abc import Sequence
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

import polars as pl
from loguru import logger


def log_and_exit(message, code=1):
    logger.error(message)
    sys.exit(code)
    
    
def filter_count(df, **kwargs):
    filters = {
        kwargs['xc']: kwargs['x_limit'],
        kwargs['yc']: kwargs['y_limit']
    }
    if 'count_NTC' in df.columns and 'ntc_limit' in kwargs:
        filters['count_NTC'] = kwargs['ntc_limit']
        
    expressions, parts = [], []
    for key, (low, high) in filters.items():
        if low is None and high is None:
            continue
        
        if low is not None:
            expressions.append(pl.col(key) >= low)
        if high is not None:
            expressions.append(pl.col(key) <= high)
        
        if low is not None and high is not None:
            parts.append(f"{low} <= {key} <= {high}")
        elif low is not None:
            parts.append(f"{key} >= {low}")
        else:
            parts.append(f"{key} <= {high}")
            
    expr = pl.all_horizontal(expressions) if expressions else None
    if expr is not None:
        logger.debug(f'Filtering {df.height:,} compounds with {len(parts):,} filters:')
        parts = parts if len(parts) == 1 else [f'({p})' for p in parts]
        fs = ' & '.join(parts)
        logger.debug(f'  {fs}')
        df = df.filter(expr)
        logger.debug(f'Retained {df.height:,} compounds passed {len(parts):,} filters')
    return df


def _enrich(row, xc, xs, ntc):
    count, score = row[xc], row[xs]
    count_ntc, score_ntc = row[f'count_{ntc}'], row[f'zscore_{ntc}']
    if score_ntc > 1 or (score_ntc > 0.5 and count_ntc > 5) or count_ntc > 20:
        return ''

    if row['nHH'] > 10:
        return ''
    if (score >= 1 and 10 <= count < 300) or (score >= 3 and 10 <= count < 100):
        return 'Weak'
    elif score > 10 and count >= 300:
        return 'Strong'
    return ''


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
        
    ntc = [name for name in names if name.lower() == 'ntc']
    ntc = ntc[0] if ntc else ''
    
    options['x'] = x
    options['y'] = y
    options['ntc'] = ntc
    options['xc'], options['xs'] = f'count_{x}', f'zscore_{x}'
    options['yc'], options['ys'] = f'count_{y}', f'zscore_{y}'
    options['x_label'] = kwargs.get('x_label', f'{x} (zscore)')
    options['y_label'] = kwargs.get('y_label', f'{y} (zscore)')
    
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
    options['x_limit'] = kwargs.get('x_limit', OPTIONS['x_limit']['default'])
    options['y_limit'] = kwargs.get('y_limit', OPTIONS['y_limit']['default'])
    options['ntc_limit'] = kwargs.get('ntc_limit', OPTIONS['ntc_limit']['default'])
    options['dpi'] = kwargs.get('dpi', OPTIONS['dpi']['default'])
    
    options['outdir'] = Path(kwargs.get('outdir', table.parent.absolute()))
    if not options['outdir'].exists():
        options['outdir'].mkdir(exist_ok=True, parents=True)
    return df, options


def process_data(df, **kwargs):
    logger.debug(f'Identifying unique compounds from {df.height:,} features ...')
    df = (df.sort(kwargs['xs'], descending=True).unique(subset=['library', 'c1_smiles', 'c2_smiles', 'c3_smiles']))
    logger.debug(f'Retained {df.height:,} unique features')
    
    df = filter_count(df, **kwargs)
    
    df = df.with_columns(
        pl.when(~pl.col('history_hits').is_null())
        .then(pl.col('history_hits').str.count_matches(';') + 1)
        .otherwise(0)
        .alias('nHH')
    )
    
    colors = {str(i): kwargs['colors'][i // 3] for i in range(7)}
    df = df.with_columns(pl.col('axis').cast(pl.String).replace(colors).alias('color'))
    
    xc, ntc = kwargs['xc'], kwargs['ntc']
    xs, ys = kwargs['xs'], kwargs['ys']
    if ntc:
        df = df.with_columns(
                pl.struct([kwargs['xc'], kwargs['xs'], f'count_{ntc}', f'zscore_{ntc}', 'history_hits', 'nHH'])
                .map_elements(lambda row: _enrich(row, xc, xs, ntc), return_dtype=pl.String)
                .alias('enrichment')
        )
    else:
        df = df.with_columns(pl.lit('').alias('enrichment'))
        
    dd = df.group_by('library').agg(
            pl.col(xs).min().alias(f'x_min'),
            pl.col(xs).max().alias(f'x_max'),
            pl.col(ys).min().alias(f'y_min'),
            pl.col(ys).max().alias(f'y_max'),
    )
    limits = {
        'min': min(dd['x_min'].min(), dd['y_min'].min()),
        'max': max(dd['x_max'].max(), dd['y_max'].max())
    }
    logger.debug(f"Dataset {limits}")
    
    for row in dd.to_dicts():
        library = row.pop('library')
        limits[library] = row
        logger.debug(f'{library} {row}')
    
    
    return df, limits


def overview_plot(df, limits, **kwargs):
    image = kwargs['outdir'].joinpath('overview.png')
    libraries = kwargs['libraries']
    num, cols = len(libraries), 6
    rows, mod = divmod(num, cols)
    rows = rows + 1 if mod else (rows or 1)

    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(18, 9))
    for i, ax in enumerate(axes.flatten()):
        row, col = divmod(i, cols)
        try:
            library = libraries[i]
        except IndexError:
            axes[row - 1][col].xaxis.set_major_locator(MaxNLocator(nbins=4))
            
            ax.set_xticks([]), ax.set_xticklabels([])
            ax.set_yticks([]), ax.set_yticklabels([])
            ax.get_xaxis().set_visible(False), ax.get_yaxis().set_visible(False)
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            continue
        dl = df.filter(pl.col('library') == library)
        enrichment = dl['enrichment'].to_list()
        if 'Weak' in enrichment:
            color = '#FF7375'
        elif 'Strong' in enrichment:
            color = '#E60003'
        else:
            color = 'black'
        
        ax.scatter(dl[kwargs['xs']], dl[kwargs['ys']], c=dl['color'],
                   s=60, edgecolors='lightgray', lw=0.3, zorder=2)
        
        ax.set_title(library, fontweight='bold', color=color, y=0.88)
        ax.axvline(x=0, color='lightgray', lw=0.5)
        ax.axhline(y=0, color='lightgray', lw=0.5)
        ax.set_xlim(limits['min'], limits['max'])
        ax.set_ylim(limits['min'], limits['max'])
        
        if row == rows - 1:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        else:
            ax.set_xticks([])
        if col == 0:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        else:
            ax.set_yticks([])
    
    fig.text(0.5, 0.01, kwargs['x_label'], ha='center', va='center', fontsize=12)
    fig.text(0.02, 0.5, kwargs['y_label'], ha='center', va='center', rotation='vertical', fontsize=12)
    
    fig.subplots_adjust(left=0.05, right=0.99, top=0.97, bottom=0.05, wspace=0.06, hspace=0.05)
    fig.savefig(image, dpi=kwargs['dpi'])
    logger.debug(f'Overview plot was saved to {image}\n')


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
    unique, limits = process_data(df, **options)
    overview_plot(unique, limits, **options)
    
    
if __name__ == '__main__':
    epv(Path('__file__').resolve().parent / 'data/fake.data.tsv.gz', verbose=True)
    