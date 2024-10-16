# Pandas Implementation (See humans2.py for polars Implementation)
import pandas as pd
from pathlib import Path


def _classify_motif(motif):
   if motif.startswith('CG'): return 'CG'
   elif motif[0] == 'C' and motif[2] == 'G': return 'CHG'
   elif motif[0] == 'C' and motif[1] != 'G' and motif[2] != 'G': return 'CHH'
   else:return 'Other'


def read_human_tsv(fn: str| Path) -> pd.DataFrame:


   df = pd.read_csv(fn, sep='\t', low_memory=False)
   return (df.loc[lambda x: x['Chrom'].str.match(r'^\d+$|^X$|^Y$', na=False)]
              .astype({'Chrom': 'category'})
              .assign(Chrom=lambda x: x['Chrom']
                      .cat.set_categories([str(i) for i in range(1, 23)] + ['X', 'Y'], ordered=True))
              .assign(Motif_type=lambda x: x['Motif'].apply(_classify_motif)))

def human_sort_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort and transform human genomic DataFrame for continuous genome-wide visualization.
    This function is specifically designed for human genomic data. It:
    1. Calculates new 'Pos' values to be continuous across all chromosomes.
    2. Adds 'Pos_chrom' column, containing original positions.
    """
    df['Pos_chrom'] = df['Pos'].copy()
    tpos_max = 0
    current_chrom = None
    new_pos = []
    for chrom, pos in zip(df["Chrom"], df["Pos_chrom"]):
        if chrom != current_chrom:
            current_chrom = chrom
            tpos_max = max(tpos_max, max(new_pos) if new_pos else 0)
        new_pos.append(tpos_max + pos)
    df["Pos"] = new_pos
    return df

def human_conv_unconv_df(df,depth=1):
    return(df[df['Depth'] >= depth]
                .groupby(['Sample'], observed=True)
                .agg({'Converted': 'sum', 'Unconverted':'sum','Depth':'sum'})
                .assign(Ratio_conv=lambda x: x['Converted'] / x['Depth'])
                .assign(Ratio_unconv=lambda x: x['Unconverted'] / x['Depth'])
                .reset_index()
                [['Sample', 'Ratio_conv','Ratio_unconv']])

def human_motif_df(df, depth=1):
    "Depth >= 2"
    return(df[df['Depth'] >= depth]
           .groupby(['Sample','Motif_type'],observed=True)
           .agg({'Ratio_conv': 'mean', 'Ratio_unconv': 'mean'})
           .reset_index())

def human_concat_dfs(in_path, samples, function_df, suffix='tsv.gz'):
    fn = fname(in_path, samples[0],suffix)
    df = read_human_tsv(fn)
    dfs = function_df(df)
    if len(samples) == 1: return dfs
    for sample in samples[1:]:
        fn = fname(in_path, sample, suffix)
        df = read_human_tsv(fn)
        df = function_df(df)
        dfs = pd.concat([dfs,df],ignore_index=True)
    return(dfs)

