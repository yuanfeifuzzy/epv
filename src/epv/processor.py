import polars as pl
from loguru import logger


def filter_count(df, **kwargs):
    filters = {}
    if any(kwargs['x_limit']) and kwargs['x_limit'][0] > 1:
        filters[kwargs['xc']] = kwargs['x_limit']
    
    if any(kwargs['y_limit']) and kwargs['y_limit'][0] > 1:
        filters[kwargs['yc']] = kwargs['y_limit']
    
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
    return df.filter(
        (pl.col('axis') == 6) & (pl.col(xs) >= kwargs['top_hits'][0]) & (pl.col(yc) <= kwargs['top_hits'][1])).sort(
            kwargs['xs'], descending=True)


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
    # df = filter_count(df, **kwargs)
    
    logger.debug(f'Identifying unique compounds from {df.height:,} features ...')
    keys = ['library', 'C1_SMILES', 'C2_SMILES', 'C3_SMILES']
    df = (
        df.with_columns(pl.len().over(keys).alias('encodings')).sort(kwargs['xs'], descending=True).unique(subset=keys))
    logger.debug(f'Retained {df.height:,} unique features')
    
    df = df.with_columns(pl.when(~pl.col('history_hits').is_null()).then(
        pl.col('history_hits').str.count_matches(';') + 1).otherwise(0).alias('nHH'))
    
    colors = {str(i): kwargs['colors'][i // 3] for i in range(7)}
    df = df.with_columns(pl.col('axis').cast(pl.String).replace(colors).alias('color'))
    
    xc, ntc = kwargs['xc'], kwargs['ntc']
    xs, ys = kwargs['xs'], kwargs['ys']
    if ntc:
        df = df.with_columns(
                pl.struct([kwargs['xc'], kwargs['xs'], f'count_{ntc}', f'zscore_{ntc}', 'nHH']).map_elements(
                    lambda row: _enrich(row, xc, xs, ntc), return_dtype=pl.String).alias('enrichment'))
    else:
        df = df.with_columns(pl.lit('').alias('enrichment'))
    
    return df