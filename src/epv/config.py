import io
import sys
from pathlib import Path
from functools import partial
from datetime import datetime
from collections.abc import Sequence

from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw, Descriptors

import polars as pl
from PIL import Image
from loguru import logger
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, DrawingArea, TextArea

OPTIONS = {
    'x': {'help': 'A column name will be assigned as x variable'},
    'y': {'help': 'A column name will be assigned as y variable'},
    'count_columns': {'help': 'Names of columns contain count data', 'type': Sequence},
    'score_columns': {'help': 'Names of columns contain zscore data', 'type': Sequence},
    'x_label': {'help': 'Label string for x-axis'},
    'y_label': {'help': 'Label string for y-axis'},
    'libraries': {'help': 'Names of libraries need to be processed', 'type': Sequence},
    'colors': {'help': 'A sequence of 3 hex color codes for mono-, di-, and tri-sython',
               'default': ('#0000ad', '#00ff33', '#ffcc00'), 'type': Sequence},
    'cards': {'help': 'Maximum number of cards to plot for each library', 'default': 6, 'type': int},
    'show_smiles': {'help': 'Show final SMILES on compound card', 'default': True, 'type': bool},
    'date': {'help': 'A date string or defaults to today'},
    'x_limit': {'help': 'Lower and upper bound of x axis variable count', 'default': [1, None], 'type': Sequence},
    'y_limit': {'help': 'Lower and upper bound of y axia variable count', 'default': [1, None], 'type': Sequence},
    'ntc_limit': {'help': 'Lower and upper bound of NTC count', 'default': [0, 20], 'type': Sequence},
    'dpi': {'help': 'The DPI of the output plot', 'default': 300, 'type': int},
    'figure_width': {'help': 'The width of figure in inch', 'default': 18, 'type': int},
    'figure_height': {'help': 'The height of figure in inch', 'default': 9, 'type': int},
    'circle_sizes': {'help': 'The circle sizes of overview and individual library scatter plots',
                     'default': [50, 500], 'type': Sequence},
    'top_hits': {'help': 'The minimum x score and maximum y count to determine top hits',
                 'default': [1.0, 0], 'type': Sequence},
    'outdir': {'help': 'Path string to a output directory'},
    'verbose': {'help': 'Enable verbose mode to emit verbose messages', 'type': bool, 'default': False}
}

COLUMNS = [
    'library',
    'axis',
    'c1',
    'c2',
    'c3',
    'c1_smiles',
    'c2_smiles',
    'c3_smiles',
]


def config_logger(level):
    logger.remove()
    formatter = '<level>[{time:YYYY-MM-DD HH:mm:ss}] {message}</level>'
    logger.add(sys.stdout, format=formatter, level=level)
    
    
def log_and_exit(message, code=1):
    logger.error(message)
    sys.exit(code)


def validate_data_and_options(table, **kwargs):
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
    
    smiles_columns = {c: c.upper() for c in columns if 'SMILES' in c.upper()}
    df = df.rename(smiles_columns)
    options['smiles_columns'] = list(smiles_columns.values())
    
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
        options['libraries'] = sorted(libraries)
    else:
        options['libraries'] = sorted(libs)
    
    colors = kwargs.get('colors', OPTIONS['colors']['default'])
    if not isinstance(colors, Sequence):
        log_and_exit(f'Invalid colors: {colors}. Colors must be a sequence of 3 hex color code strings')
    colors = list(colors)
    if len(colors) != 3:
        log_and_exit(f'Invalid colors: {colors}. Colors must be a sequence of 3 hex color code strings')
    options['colors'] = colors
    
    options['cards'] = kwargs.get('cards', OPTIONS['cards']['default'])
    options['show_smiles'] = kwargs.get('show_smiles', OPTIONS['show_smiles']['default'])
    options['date'] = kwargs.get('date', datetime.now().strftime('%m%d%Y'))
    options['x_limit'] = kwargs.get('x_limit', OPTIONS['x_limit']['default'])
    options['y_limit'] = kwargs.get('y_limit', OPTIONS['y_limit']['default'])
    options['ntc_limit'] = kwargs.get('ntc_limit', OPTIONS['ntc_limit']['default'])
    options['dpi'] = kwargs.get('dpi', OPTIONS['dpi']['default'])
    options['figure_width'] = kwargs.get('figure_width', OPTIONS['figure_width']['default'])
    options['figure_height'] = kwargs.get('figure_height', OPTIONS['figure_height']['default'])
    options['circle_sizes'] = kwargs.get('circle_sizes', OPTIONS['circle_sizes']['default'])
    
    xs, yc = kwargs.get('top_hits', OPTIONS['top_hits']['default'])
    try:
        options['top_hits'] = [float(xs), int(yc)]
    except ValueError as e:
        log_and_exit(f'Invalid thresholds for determine top hits: {[xs, yc]}')
    
    options['outdir'] = Path(kwargs.get('outdir', table.parent.absolute()))
    if not options['outdir'].exists():
        options['outdir'].mkdir(exist_ok=True, parents=True)
    return df, options
