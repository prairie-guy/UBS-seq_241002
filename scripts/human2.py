### Polars Implemantation of human.py

import polars as pl
import pandas as pd
from pathlib import Path
from fnames import *

def pl2pd(polars_df: pl.DataFrame, chunk_size: int = 10_000_000) -> pd.DataFrame | None:
    """
    Convert a Polars DataFrame to a pandas DataFrame in chunks.
    Returns: The converted pandas DataFrame, or None if conversion fails.
    Useful for large datatrames
    """
    total_rows = polars_df.shape[0]
    pandas_dfs = []
    try:
        for i in range(0, total_rows, chunk_size):
            chunk = polars_df.slice(i, chunk_size)
            pandas_chunk = pd.DataFrame(chunk.to_dicts())
            pandas_dfs.append(pandas_chunk)

        return pd.concat(pandas_dfs, ignore_index=True)
    except Exception as e:
        print(f"Error during conversion: {e}")
        return None


# Includes Categorical Data for Chromosomes Sorting
def read_human_tsv(fn: str| Path) -> pl.DataFrame:
    def classify_motif(motif: str) -> str:
        if motif.startswith('CG'): return 'CG'
        elif motif.startswith('C') and len(motif) > 2 and motif[2] == 'G': return 'CHG'
        elif (motif.startswith('C') and len(motif) > 2 and motif[1] != 'G' and motif[2] != 'G'): return 'CHH'
        else:return 'Other'
    chrom_categories = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    return (pl.read_csv(fn, separator='\t',
            schema_overrides={
                'Chrom': pl.Utf8,
                'Ratio_conv': pl.Float64,
                'Ratio_unconv': pl.Float64})
        .filter(pl.col('Chrom').is_in(chrom_categories))
        .with_columns(
            pl.col('Chrom').cast(pl.Enum(chrom_categories)),
             Motif_type=pl.col('Motif').map_elements(classify_motif, return_dtype=pl.Utf8)))


# Excludes Categorical Data for Chromosomes Sorting
# def read_human_tsv(fn: str| Path) -> pl.DataFrame:
#     def classify_motif(motif: str) -> str:
#         if motif.startswith('CG'): return 'CG'
#         elif motif.startswith('C') and len(motif) > 2 and motif[2] == 'G': return 'CHG'
#         elif (motif.startswith('C') and len(motif) > 2 and motif[1] != 'G' and motif[2] != 'G'): return 'CHH'
#         else:return 'Other'

#     chrom_categories = [str(i) for i in range(1, 23)] + ['X', 'Y']
#     return (pl.read_csv(fn, separator='\t',
#             schema_overrides={
#                 'Chrom': pl.Utf8,
#                 'Ratio_conv': pl.Float64,
#                 'Ratio_unconv': pl.Float64})
#         .filter(pl.col('Chrom').is_in(chrom_categories))
#         .with_columns(Motif_type=pl.col('Motif').map_elements(classify_motif, return_dtype=pl.Utf8)))

def human_conv_unconv_df(df: pl.DataFrame, depth: int = 1) -> pl.DataFrame:
    return (df.filter(pl.col('Depth') >= depth)
              .group_by('Sample',maintain_order=True)
              .agg([pl.sum('Converted').alias('Converted'),
                    pl.sum('Unconverted').alias('Unconverted'),
                    pl.sum('Depth').alias('Depth')])
              .with_columns([
                  (pl.col('Converted') / pl.col('Depth')).alias('Ratio_conv'),
                  (pl.col('Unconverted') / pl.col('Depth')).alias('Ratio_unconv')])
              .select(['Sample', 'Ratio_conv', 'Ratio_unconv']))


def human_motif_df(df: pl.DataFrame, depth: int = 1) -> pl.DataFrame:
    return (df.filter(pl.col('Depth') >= depth)
              .group_by(['Sample', 'Motif_type'],maintain_order=True)
              .agg([pl.sum('Converted').alias('Converted'),
                    pl.sum('Unconverted').alias('Unconverted'),
                    pl.sum('Depth').alias('Depth')])
              .with_columns([
                  (pl.col('Converted') / pl.col('Depth')).alias('Ratio_conv'),
                  (pl.col('Unconverted') / pl.col('Depth')).alias('Ratio_unconv')])
              .select(['Sample','Depth', 'Motif_type', 'Ratio_conv', 'Ratio_unconv']))


def human_sort_df(df: pl.DataFrame) -> pl.DataFrame:
    df = df.with_columns(pl.col('Pos').alias('Pos_chrom'))
    tpos_max, current_chrom, new_pos = 0, None, []
    for chrom, pos in zip(df["Chrom"].to_list(), df["Pos_chrom"].to_list()):
        if chrom != current_chrom:
            current_chrom = chrom
            tpos_max = max(tpos_max, max(new_pos) if new_pos else 0)
        new_pos.append(tpos_max + pos)
    return df.with_columns(pl.Series(new_pos).alias('Pos'))


def human_concat_dfs(in_path, samples, function_df, suffix='pq') -> pl.DataFrame:
    fn = fname(in_path, samples[0],suffix)
    df = read_human_tsv(fn) if suffix == 'tsv.gz' else pl.read_parquet(fn)
    dfs = function_df(df)
    if len(samples) == 1: return dfs
    for sample in samples[1:]:
        fn = fname(in_path, sample, suffix)
        df = read_human_tsv(fn) if suffix == 'tsv.gz' else pl.read_parquet(fn)
        df = function_df(df)
        dfs = pl.concat([dfs, df])
    return dfs

