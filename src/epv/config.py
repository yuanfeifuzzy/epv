from collections.abc import Sequence
from argparse import BooleanOptionalAction

OPTIONS = {
    'x': {'help': 'A column name will be assigned as x variable'},
    'y': {'help': 'A column name will be assigned as y variable'},
    'count_columns': {'help': 'Names of columns contain count data', 'type': Sequence},
    'score_columns': {'help': 'Names of columns contain zscore data', 'type': Sequence},
    'x_label': {'help': 'Label string for x-axis'},
    'y_label': {'help': 'Label string for y-axis'},
    'libraries': {'help': 'Names of libraries need to be processed', 'type': Sequence},
    'colors': {'help': 'A sequence of 3 hex color codes for mono-, di-, and tri-sython',
               'default': ('#0000ad', '#00ff33', '#ffff00'), 'type': Sequence},
    'cards': {'help': 'Maximum number of cards to plot for each library', 'default': 6, type: int},
    'date': {'help': 'A date string or defaults to today'},
    'outdir': {'help': 'Path string to a output directory'},
    'verbose': {'help': 'Enable verbose mode to emit verbose messages', 'type': BooleanOptionalAction, 'default': False}
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