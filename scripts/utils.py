#!/usr/bin/env python3

# `utils.py`
#
# C. Bryan Daniels, cdaniels@nandor.net
# 03/17/2024

"""
Utilities to be used with Jupyter
"""

from find_cmers import *
from fnames import *
from reference import *

from pathlib import Path
from typing import Dict, Optional, Tuple, List, Union
import os, sys, re, subprocess, itertools, colorsys, importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def reload_module(module_name: str) -> None:
    """Reloads the specified module and updates the caller's global namespace."""
    module = importlib.import_module(module_name)
    importlib.reload(module)
    caller_globals = sys._getframe(1).f_globals
    caller_globals.update({name: getattr(module, name) for name in dir(module) if not name.startswith('_')})

def adaptive_formatter(threshold=1e10):
    'Set threshold higher than 1e9 to handle G, M, K formatting for matplotlib plt'
    def format_number(x, pos):
        if abs(x) >= threshold:
            return f'{x:.1e}'
        elif abs(x) >= 1e9:
            return f'{x/1e9:.1f}G'
        elif abs(x) >= 1e6:
            return f'{x/1e6:.1f}M'
        elif abs(x) >= 1e3:
            return f'{x/1e3:.0f}K'
        else:
            return f'{x:.0f}'
    return plt.FuncFormatter(format_number)


# Windows Functions for setting windows
def windows_unconverted(df, size=100):
    "Windows for UBS-seq"
    return (
        df.assign(Pos=lambda x: x.Pos // size * size)
        .groupby(["Chrom", "Pos"])
        .agg({"Unconverted": "sum", "Converted": "sum"})
        .assign(Ratio=lambda x: x.Unconverted / (x.Converted + x.Unconverted))
        .reset_index())

def windows_converted(df, size=100):
    "Windows for BAT-seq"
    return (
        df.assign(Pos=lambda x: x.Pos // size * size)
        .groupby(["Chrom", "Pos"])
        .agg({"Unconverted": "sum", "Converted": "sum"})
        .assign(Ratio=lambda x: x.Converted / (x.Converted + x.Unconverted))
        .reset_index())

def windows_depth(df, size=100):
    "Windows for Depth Coverage"
    return (
        df.assign(Window=lambda x: x.Pos // size * size)
        .groupby(["Strand", "Window"])
        .agg({"Depth": "mean", "Pos": "first"})  # Use mean depth and first position in window
        .reset_index()
        .drop(columns=["Window"]))

def reorder_df(df: pd.DataFrame, sample: list, key: str = 'SampleID') -> pd.DataFrame:
    """Reorder DataFrame rows to match the order of values in the given list.

    Args:
        df: Input DataFrame
        sample: List of values to order by
        key: Column name to use for ordering (default: 'SampleID')

    Returns:
        Reordered DataFrame
    """
    return df.set_index(key).loc[sample].reset_index()

def mk_read_counts(out_paths,samples):
    """
    Takes out_path::DictList or instance of Path. Only expects (file_R1.fq.gz,file_R2.fq.gz) or file.bam
    Creates a out_path/'prefix'_read_counts.gz file
    """
    nc = os.cpu_count()
    if isinstance(out_paths,Path): out_paths = [out_paths]
    (tail,ext) = ('','bam') if list(out_paths[0].glob('*.bam')) else ('_R1','fq.gz')
    for out_path in out_paths:
        prefix = str(out_path).split('_')[0]
        with open(fname(out_path,f'{prefix}_read_counts', 'tsv'),'w') as fh:
            print(f'SampleID\tCount',file=fh)
            for sample in samples:
                if tail == '_R1': cmd = ["samtools", "view", "-@", str(nc), "-c", fname(out_path, sample + tail, ext)]
                else:             cmd = ["samtools", "view","-@", str(nc), "-F4","-F16","-F256", "-c", fname(out_path, sample + tail, ext)]
                result = subprocess.run(cmd, capture_output=True, text=True)
                c = int(result.stdout.strip())
                print(f'{sample}\t{c}',file=fh)

def get_read_counts(out_paths):
    """
    Takes an out_paths::fname.ListDict or a Path instance and samples
    It looks for a file of the form 'some_prefix'_read_counts.tsv with headers: SampleID, Counts
    Warning: It will fail if the title does not have the suffix or the headers differ
    Returns a merged df with SampleID, Count(ref[0]), Count(ref[1])
    Example: get_read_counts(out_paths, 'post_dedup_filter_read_counts.tsv')
    """
    nc = os.cpu_count()
    if isinstance(out_paths, Path): out_paths = [out_paths]
    try:
        return (pd.concat([
            pd.read_csv(list(out_path.glob('*_read_counts.tsv'))[-1], sep='\t')
            .rename(columns={'Count': str(out_path).split('_')[-1]})
            for out_path in out_paths])
            .groupby('SampleID')
            .first()
            .reset_index()
            .apply(lambda col: col.astype(int) if col.name != 'SampleID' else col))
    except:
        print("error: Files missing of form 'some_prefix'_read_counts.tsv with headers: SampleID, Counts")
        return None

def sanitize(title: str) -> str:
    """Sanitize string by replacing specific chars with underscores."""
    sanitized = re.sub(r'[^a-zA-Z0-9\s]', '_', title)
    sanitized = re.sub(r'\s+', '_', sanitized)
    sanitized = re.sub(r'_+', '_', sanitized)
    return sanitized.strip('_')


def generate_bar_graph(
    df: pd.DataFrame,
    value_col: str,
    label_col: str = 'SampleID',
    groups: Optional[Dict[str, List[int]]] = None,
    title: Optional[str] = None,
    ylabel: Optional[str] = None,
    yrange: Union[float, Tuple[float, float]] = 100,
    size: Tuple[int, int] = (22, 6),
    digits: int =1
) -> None:
    """
    Generate a bar graph from the given dataframe.

    Parameters:
        df (pd.DataFrame): Input dataframe

        value_col (str): Column name for y-axis values

        label_col (str): Column name for x-axis labels (default: 'SampleID')

        groups (Optional[Dict[str, List[int]]]): Dictionary of group names and their member row indices (0-indexed) (default: None)

        title (Optional[str]): Title of the graph (default: None)

        ylabel (Optional[str]): Label for y-axis (default: None, uses value_col name if not provided)

        yrange (Union[float, Tuple[float, float]]): Maximum value or range for y-axis (default: 100)

        size (Tuple[int, int]): Figure size (default: (22, 6))

        digits (int): Truncated label above bar to this number of digits

    Returns:
        None (displays and saves the graph)

    Note:
        - If a 'std' column is present in the DataFrame, it will be used to display error bars.
    """
    if label_col not in df.columns:
        raise ValueError(f"Label column '{label_col}' not found in DataFrame.")
    if value_col not in df.columns:
        raise ValueError(f"Value column '{value_col}' not found in DataFrame.")
    if df[label_col].isnull().any():
        raise ValueError(f"Column '{label_col}' contains missing values.")
    if df[value_col].isnull().any():
        raise ValueError(f"Column '{value_col}' contains missing values.")

    plt.figure(figsize=size)
    if groups:
        groups = dict(groups)
        base_palette = sns.color_palette("husl", len(groups))
        group_colors = {name: base_palette[i] for i, name in enumerate(groups)}
        df_plot = pd.DataFrame()
        for group_name, indices in groups.items():
            group_data = df.iloc[indices].copy()
            group_data['group'] = group_name
            df_plot = pd.concat([df_plot, group_data])
        df_plot = df_plot.reset_index(drop=True)
    else:
        df_plot = df.copy()
        df_plot['group'] = ''
        group_colors = None

    ax = sns.barplot(
        x=label_col,
        y=value_col,
        hue='group',
        data=df_plot,
        palette=group_colors,
        dodge=False
    )
    if 'std' in df_plot.columns:
        ax.errorbar(
            x=range(len(df_plot)),
            y=df_plot[value_col],
            yerr=df_plot['std'],
            fmt='none',
            c='black',
            capsize=5
        )
    plt.title(title if title else '', fontsize=20)
    plt.xlabel('')
    plt.ylabel(ylabel if ylabel else value_col)
    if isinstance(yrange, tuple):
        plt.ylim(yrange[0], yrange[1])
    else:
        plt.ylim(0, yrange)
    plt.xticks(rotation=45, ha='right')
    for i, v in enumerate(df_plot[value_col]):
        label = f'{v:.{digits}f}'
        ax.text(i, v, label, ha='center', va='bottom')
    if groups:
        ax.legend(title='Groups', loc='upper right')
    else:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
    plt.tight_layout()
    title = sanitize(title)
    save_path = Path('figures') / f'{title}.png'
    plt.savefig(save_path, format='png', bbox_inches='tight')
    plt.show()


def agg_replicates(df: pd.DataFrame,
                  group_size_or_list: Union[int, List[List[int]]],
                  sample_col: str = 'SampleID',
                  value_col: str = 'ac') -> pd.DataFrame:
    """
    Aggregate groups of rows in a DataFrame, computing mean and standard deviation.

    Parameters:
    df (pd.DataFrame): Input DataFrame
    group_size_or_list (int or list of lists): If int, divides rows into groups of that size.
                                               If list of lists, groups rows according to the lists.
    sample_col (str): Column name for sample identifiers (default: 'SampleID')
    value_col (str): Column name for value_col to be averaged (default: 'ac')

    Returns:
    pd.DataFrame: Aggregated DataFrame with columns 'SampleID', 'mean', and 'std'
    """
    if (not isinstance(df, pd.DataFrame) or not isinstance(group_size_or_list, (int, list))
        or (isinstance(group_size_or_list, int) and len(df) % group_size_or_list != 0)):
        return "agg_replicates(df: pd.DataFrame, group_size_or_list: Union[int, List[List[int]]], \
                sample_col: str = 'SampleID', value_col: str = 'ac') -> pd.DataFrame"

    df = df.assign(**{value_col: pd.to_numeric(df[value_col], errors='coerce')})
    grouping = (np.arange(len(df)) // group_size_or_list if isinstance(group_size_or_list, int)
                else [i for i, group in enumerate(group_size_or_list) for _ in group])
    return (df.groupby(grouping)
              .agg({sample_col: 'first', value_col: ['mean', 'std']})
              .reset_index(drop=True)
              .pipe(lambda x: pd.DataFrame({
                  'SampleID': x[sample_col]['first'],
                  'mean': x[value_col]['mean'],
                  'std': x[value_col]['std']})))
