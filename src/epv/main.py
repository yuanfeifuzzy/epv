import functools
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import polars as pl
from drawer import Card, position_cards
from config import OPTIONS, logger, config_logger, log_and_exit, validate_data_and_options
    
    
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


def find_top_hits(df, **kwargs):
    xs, yc = kwargs['xs'], kwargs['yc']
    return df.filter((pl.col('axis') == 6) &
                     (pl.col(xs) >= kwargs['top_hits'][0]) &
                     (pl.col(yc) <= kwargs['top_hits'][1])).sort(kwargs['xs'], descending=True)


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


def process_data(df, **kwargs):
    logger.debug(f'Identifying unique compounds from {df.height:,} features ...')
    df = (df.sort(kwargs['xs'], descending=True).unique(subset=['library', 'C1_SMILES', 'C2_SMILES', 'C3_SMILES']))
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

    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(kwargs['figure_width'], kwargs['figure_height']))
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
                   s=kwargs['circle_sizes'][0], edgecolors='#6f6f6f', lw=0.3, zorder=2)
        
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


def adjust_limits_for_marker(ax, limit: dict, s=500, offset_pts=10):
    # 1. Calculate marker radius in points
    # Area s = pi * r^2  => r = sqrt(s/pi)
    radius_pts = np.sqrt(s / np.pi)
    
    # 2. Total physical padding needed (radius + your offset)
    total_pad_pts = radius_pts + offset_pts
    
    # 3. Convert points to inches (72 points = 1 inch)
    pad_inches = total_pad_pts / 72.0
    
    # 4. Get the current axes size in inches
    # This requires the renderer to have run, or we use the figure size
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    ax_width_in, ax_height_in = bbox.width, bbox.height
    
    # 5. Calculate data ranges
    x_min, x_max, y_min, y_max = limit['x_min'], limit['x_max'], limit['y_min'], limit['y_max']
    x_range, y_range = x_max - x_min, y_max - y_min
    
    # 6. Convert physical padding to data units
    # (Padding_Inches / Total_Inches) * Total_Data_Range
    x_pad_data = (pad_inches / ax_width_in) * x_range
    y_pad_data = (pad_inches / ax_height_in) * y_range
    
    ax.set_xlim(x_min - x_pad_data, x_max + x_pad_data)
    ax.set_ylim(y_min - y_pad_data, y_max + y_pad_data)
    
    
def set_title_legend(ax, library, **kwargs):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="8%", pad=0)
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.add_artist(AnchoredText(library, loc='upper right', prop=dict(size=14, fontweight='bold'), frameon=False))
    
    c1, c2, c3 = kwargs['colors']
    options = {'marker': 'o', 'markersize': 15, 'markeredgecolor': '#6f6f6f', 'color': 'w'}
    legends = [Line2D([0], [0], label='Mono-', markerfacecolor=c1, **options),
               Line2D([0], [0], label='Di-', markerfacecolor=c2, **options),
               Line2D([0], [0], label='Tri-Sython', markerfacecolor=c3, **options)]
    cax.legend(handles=legends, loc='center left', bbox_to_anchor=(0, 0.5), frameon=False, ncol=3, fontsize=10)
    
    
def compound_cards(df, du, ax, fig, **kwargs):
    n, tops = kwargs['cards'], find_top_hits(du, **kwargs)
    if n and not tops.is_empty():
        tops, cards = tops.slice(0, n).reverse(), []
        xs, width = position_cards(ax, tops[kwargs['xs']].to_list())
        y = du[kwargs['ys']].max()
        
        for x, row in zip(xs, tops.iter_rows(named=True)):
            smiles_columns, show_smiles = kwargs['smiles_columns'], kwargs['show_smiles']
            smiles_columns = [s for s in smiles_columns if s != 'SMILES']
            smiles = ['SMILES'] + smiles_columns if show_smiles else smiles_columns
            smiles = [row[s] for s in smiles]
            data = {name: f'{row[cc]} ({row[sc]})'
                    for name, cc, sc in zip(kwargs['names'], kwargs['count_columns'], kwargs['score_columns'])}
            data['nHH'] = row['nHH']
            point = (row[kwargs['xs']], row[kwargs['ys']])
            card = Card(row.get('compound', ''), smiles, data, point, x, y, width, ax, fig)
            card.render()
    else:
        if n:
            logger.debug('No top hits was found, no compound card will be drawn')
        else:
            logger.debug('No top hits was found, no compound card will be drawn')
    
    
def library_plot(df, du, library, limits, **kwargs):
    logger.debug(f'Plotting library {library} ...')
    dd = du.filter(pl.col('library') == library)
    if dd.height == 0:
        logger.debug(f'No compounds passed filters, skip plotting library {library}')
    else:
        limit = limits[library]
        logger.debug(f'Setting figure limits to {limit}')
        ax = kwargs.get('ax', None)
        fig, ax = (plt, ax) if ax else plt.subplots(figsize=(kwargs['figure_width'], kwargs['figure_height']))
        ax.scatter(dd[kwargs['xs']], dd[kwargs['ys']], c=dd['color'], s=kwargs['circle_sizes'][1], edgecolors='#6f6f6f',
                   lw=0.25, zorder=2, clip_on=False)
            
        ax.set_xlabel(kwargs['x_label'], fontsize=12), ax.set_ylabel(kwargs['y_label'], fontsize=12)
        adjust_limits_for_marker(ax, limit, s=kwargs['circle_sizes'][1])
        set_title_legend(ax, library, **kwargs)
        compound_cards(df, dd, ax, fig, **kwargs)
        
        image = kwargs['outdir'].joinpath(f'{library}.png')
        fig.savefig(image, dpi=kwargs['dpi'])
        logger.debug(f'Library plot was saved to {image.name}\n')
        plt.close()


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
    unique, limits = process_data(df, **options)
    # overview_plot(unique, limits, **options)
    for library in options['libraries']:
        library_plot(df, unique, library, limits, **options)
    
    
if __name__ == '__main__':
    epv(Path('__file__').resolve().parent / 'data/fake.data.tsv.gz', verbose=True, libraries=['qDOS18_1'],
        top_hits=[0.004, 0])
    # position_cards(0, 10, [8, 9])
    