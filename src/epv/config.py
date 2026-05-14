from collections.abc import Sequence

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
    'date': {'help': 'A date string or defaults to today'},
    'x_limit': {'help': 'Lower and upper bound of x axis variable count', 'default': [1, None], 'type': Sequence},
    'y_limit': {'help': 'Lower and upper bound of y axia variable count', 'default': [1, None], 'type': Sequence},
    'ntc_limit': {'help': 'Lower and upper bound of NTC count', 'default': [None, 20], 'type': Sequence},
    'dpi': {'help': 'The DPI of the output plot', 'default': 300, 'type': int},
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