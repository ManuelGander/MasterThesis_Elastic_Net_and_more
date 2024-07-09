import os
import re
import sys
import argparse
from pathlib import Path
import warnings
import logging
import math
from itertools import product
from typing import Dict, List

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from ppca import PPCA
from umap import UMAP

from . import utils
from . import config
from . import basket_scoring
from . import check_results_quality

pd.set_option('display.max_rows', None)

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)

ALL_COLORS = ['silver', 'dodgerblue', 'darkorange', 'limegreen', 'mediumorchid', 'sienna', 'darkturquoise', 'darkkhaki',
              'blueviolet', 'seagreen', 'orange', 'tan', 'slategrey', 'olive', 'lightpink', 'darkslategray',
              'mediumvioletred', 'khaki', 'powderblue', 'lightsalmon', 'olivedrab', 'firebrick',
              'lawngreen', 'steelblue', 'indigo', 'linen', 'springgreen', 'gold', 'darkred', 'lightgreen', 'pink',
              'yellowgreen']
ALL_SYMBOLS = ['o', 'v', '^', 's', '+', 'd', 'x', 'p', 'P', '<', '>']


# MT: rename function
def do_pca(selected_proteins: List[str],
           results_folder: str,
           plot_types: List[str],
           sample_annot: str,
           metadata_annot: str,
           data_type: str,
           min_sample_occurrence_ratio: float = 0.5,
           umap: bool = False,
           include_reference_channels: bool = False,
           include_replicates: bool = False,
           plot: bool = False):
    """
    data_type = fp or pp
    """
    # load sample and metadata
    sample_annot_df = pd.read_csv(sample_annot)
    sample_annot_df.columns = sample_annot_df.columns.str.strip()
    # sample_annot_df = sample_annot_df[sample_annot_df['QC'] == 'passed']

    meta_annot_df = pd.read_excel(metadata_annot)
    meta_annot_df = whitespace_remover(meta_annot_df)

    if include_reference_channels:
        sample_annot_df = create_ref_sample_annot(results_folder, sample_annot_df)

    # check error
    metadata_df = merge_sample_and_metadata_annots(sample_annot_df,
                                                   meta_annot_df,
                                                   remove_qc_failed=True,
                                                   keep_reference=include_reference_channels,
                                                   keep_replicates=include_replicates)

    all_principal_dfs = []
    all_principal_variances = []
    # for data_type in data_types:
    for plot_type in plot_types:  # MT: plot types are actually different data inputs, change name
        logger.info('plot type: ', plot_type)
        df = load_pca_data(results_folder, metadata_df['Sample name'], data_type, plot_type, include_reference_channels)
        if len(selected_proteins) > 2:
            final_selected = [f for f in selected_proteins if f in df.index]
            if len(final_selected) > 2:
                df = df[df.index.isin(final_selected)]

        principal_df, pca = metadata_pca(df,
                                         metadata_df,
                                         plot_type,
                                         umap,
                                         min_sample_occurrence_ratio)
        feature_importance_pc = check_results_quality.investigate_pca_feature_importance(df, pca)
        print(feature_importance_pc.head(15))
        all_principal_dfs.append(principal_df)
        if umap:
            all_principal_variances.append([])
        else:
            all_principal_variances.append(pca.var_exp)

        if plot:
            plot_folder = os.path.join(results_folder, 'Results_investigation_plots', f'{plot_type}')
            os.makedirs(plot_folder, mode=0o700, exist_ok=True)
            for metadata_column in metadata_df.columns:
                figure_name, figure_title = get_figure_name_and_title(
                    metadata_column,
                    plot_folder,
                    plot_type,
                    data_type,
                    umap,
                    include_reference_channels,
                    include_replicates,
                    min_sample_occurrence_ratio)
                plot_pca(principal_df,
                         pca,
                         metadata_df[metadata_column],
                         figure_name,
                         figure_title,
                         umap)
    return all_principal_dfs, all_principal_variances, metadata_df


def get_figure_name_and_title(metadata_column: str,
                              plot_folder: str,
                              plot_type: str,
                              data_type: str,
                              umap: bool,
                              include_reference_channels: bool,
                              include_replicates: bool,
                              min_sample_occurrence_ratio: float):
    pca_type = 'PCA'
    if umap:
        pca_type = 'UMAP'

    extra = '_'
    if include_replicates:
        extra += 'rep_'
    if include_reference_channels:
        extra += 'ref'
    figure_name = os.path.join(plot_folder,
                               f'{plot_type}_{pca_type}_{min_sample_occurrence_ratio}_{data_type}_{metadata_column}{extra}.png')
    figure_title = f'{pca_type}_{data_type}_{metadata_column}'

    return figure_name, figure_title


def load_pca_data(results_folder, samples, data_type: str, plot_type: str = 'Basket', include_reference_channels: bool = False):
    if plot_type == 'Basket':
        df = basket_scoring.read_basket_scores(results_folder, data_type)
    if plot_type == 'Kinase':
        df = read_kinase_score_file(results_folder)
    if plot_type in ['Intensity', 'Intensity_subset']:  # MT: move into a function
        file = f'annot_{data_type}.csv'
        if include_reference_channels:
            file = f'annot_{data_type}_with_ref.csv'

        index_col = utils.get_index_cols(data_type)
        df = pd.read_csv(os.path.join(results_folder, file), index_col=index_col)

        df.columns = df.columns.str.strip()
        if plot_type == 'Intensity_subset':
            # Subset remove where both basket and rtk is empty
            if 'rtk' in df.columns:
                df = df.dropna(subset=['basket', 'rtk'], how='all')
            elif 'sub_basket' in df.columns:
                df = df.dropna(subset=['basket', 'sub_basket'], how='all')
            else:
                logger.info('Error. Wrong basket scoring file. cannot find baskets of type basket & rtk or basket & sub_basket')
        df = utils.keep_only_sample_columns(df)

    # prepare data
    if plot_type in ['Basket', 'Kinase']:  # MT: include transpose as a parameter in the loading function
        df = df.transpose()

    if not include_reference_channels:  # but then also replicates are kept
        df = df.loc[:, df.columns.isin(samples)]  # MT: not sure what this does

    return df


# MT: split data loading from performing PCA
def metadata_pca(df: pd.DataFrame,
                 metadata_df: pd.DataFrame,
                 plot_type: str = 'Basket',
                 umap: bool = False,
                 min_sample_occurrence_ratio: float = 0.5):
    normalized_data, transposed_df = filter_and_normalize(df, min_sample_occurrence_ratio)

    # calculate ppca
    pca = calculate_pca(normalized_data, plot_type, umap)

    # get principal df
    principal_df = get_pca_df(pca, transposed_df, transposed_df.index, plot_type, umap)

    # merge principal components with metadata
    principal_df = principal_df.merge(metadata_df, left_on='Sample', right_on='Sample name')

    return principal_df, pca


# MT: split data loading from performing PCA
# def metadata_pca(df: pd.DataFrame,
#                  metadata_df: pd.DataFrame,
#                  plot_type: str = 'Basket',
#                  umap: bool = False,
#                  min_sample_occurrence_ratio: float = 0.5):
#     normalized_data, transposed_df = filter_and_normalize(df, min_sample_occurrence_ratio)
#
#     # calculate ppca
#     pca = calculate_pca(normalized_data, plot_type, umap)
#
#     # get principal df
#     principal_df = get_pca_df(pca, transposed_df, transposed_df.index, plot_type, umap)
#
#     # merge principal components with metadata
#     principal_df = principal_df.merge(metadata_df, left_on='Sample', right_on='Sample name', validate='one_to_one')
#
#     return principal_df, pca


# MT: move each of the reading functions to the corresponding file, here e.g. TOPAS_kinase_scoring.py
def read_kinase_score_file(results_folder):
    df = pd.read_csv(os.path.join(results_folder, 'kinase_results', f'kinase_scores.tsv'), sep='\t', index_col='PSP Kinases')
    df = df.transpose()
    df = df.drop('No. of total targets')
    return df


def get_pca_df(pca, df, samples, plot_type: str = 'Basket', umap: bool = False, target_dict: Dict = None,
               principal_df: pd.DataFrame = None, target_name: str = None):
    if not umap:
        if plot_type in ['Intensity', 'Intensity_subset']:  # MT: don't understand what this does
            pca = pca.transform()
        else:
            pca = pca.principal_components
    if target_dict is None:
        principal_df = pd.DataFrame(data=pca[:, [0, 1]], columns=['Principal component 1', 'Principal component 2'])
        principal_df = pd.concat([principal_df, pd.Series(samples)], axis=1)
        principal_df.columns = ['Principal component 1', 'Principal component 2', 'Sample']
    else:
        # targets = get_pca_targets(df, plot_type, target_dict)
        targets = get_pca_targets(df, plot_type, target_dict)
        new_columns = principal_df.columns.append(pd.Index([f'{target_name}']))
        principal_df = pd.concat([principal_df, pd.Series(targets)], axis=1)
        principal_df.columns = new_columns
    return principal_df


# MT: rename function to something more descriptive
def filter_and_normalize(df, min_sample_occurrence_ratio: float = 0.5):
    df = df.transpose()
    x = df.loc[:, df.count(axis=0) >= df.shape[0] * min_sample_occurrence_ratio]
    x = StandardScaler().fit_transform(x)  # MT: turn into dataframe?
    return x, df


# MT: rename function, it's confusing that this function can do an umap
def calculate_pca(x, plot_type, umap: bool = False):
    if umap:
        pca = do_umap(x)
    else:
        if plot_type in ['Basket', 'Kinase']: # MT: check for presence of missing values rather than check for plot_type
            pca = PCA(n_components=2)
            pca.principal_components = pca.fit_transform(x)
            pca.var_exp = pca.explained_variance_ratio_
        else:
            pca = PPCA()
            pca.fit(data=x, d=2, verbose=False)
            pca.var_exp[1] = pca.var_exp[1] - pca.var_exp[0]
    return pca


# MT: rename function
def do_umap(x, n_comp: int = 30, n_neigh: int = 5, n_epochs: int = 1000, metric: str = 'euclidean'):
    # UMAP needs complete data (no missing values) so first PPCA is used for imputation
    ppca = PPCA()
    ppca.fit(data=x, d=2, verbose=False)
    imputed_data = pd.DataFrame(ppca.data)
    reducer = UMAP(n_components=n_comp, n_epochs=n_epochs, low_memory=True, n_neighbors=n_neigh, metric=metric)
    x_trans = reducer.fit_transform(imputed_data)
    return x_trans


# MT: rename function to get_pca_colors_and_symbols
def get_pca_colors_symbols(targets):
    marker = False  # MT: no need for this variable
    if len(targets.unique()) > 8:
        marker = True
        num_combinations = math.ceil(np.sqrt(len(targets.unique())))
    else:
        num_combinations = len(targets.unique())
    colors = ALL_COLORS[0:num_combinations]
    symbols = ALL_SYMBOLS[0]
    if marker:
        symbols = ALL_SYMBOLS[0:num_combinations]
        color_symbol_tuples = []
        for shape, color in product(colors, symbols):
            color_symbol_tuples.append((shape, color))
        color_symbol_tuples = color_symbol_tuples[0:len(targets.unique())]
    else:
        # todo: make colors to color marker tuple
        color_symbol_tuples = [(color, symbols) for color in colors]
    return color_symbol_tuples


# MT: what does this do?
def get_pca_targets(df: pd.DataFrame, plot_type,
                    sample_dict=None) -> List:
    # TODO: clean up
    """Requires that fp and pp is the same """
    targets = []
    for v in df.index.tolist(): # MT: do not use v as variable name
        if v in sample_dict.keys():
            targets.append(sample_dict[v])
        # replicates
        elif '-R1' in v or '-R2' in v or '-R3' in v or '-R4' in v: # MT: use regex in case R5, R6 come up later
            temp_v = '-'.join(v.split('-')[:-1])  # MT: turn into a function
            if temp_v in sample_dict.keys():
                targets.append(sample_dict[temp_v])
        else:
            continue
    return targets


def plot_pca(principal_df, pca, targets, figure_name, figure_title, umap: bool = False):
    fig, axes = plt.subplots(1, 1, figsize=(8.27, 8.27))
    # print(targets)
    colors_symbols = get_pca_colors_symbols(targets)
    for target, marker_tuple in zip(targets.unique(), colors_symbols):
        temp_df = principal_df[principal_df[f'{targets.name}'] == target]
        plt.scatter(temp_df[['Principal component 1']],
                    temp_df[['Principal component 2']], s=30, alpha=1, c=marker_tuple[0], marker=marker_tuple[1])
    legend_size = 8
    if len(targets.unique()) > 30:
        legend_size = 6
    plt.legend(targets.unique(), prop={'size': legend_size}, loc='center right', bbox_to_anchor=(1.15, 0.5))
    if not umap:
        axes.set_xlabel('PC1 (%.0f%%)' % (np.round(pca.var_exp[0] * 100)), fontsize=19)
        axes.set_ylabel('PC2 (%.0f%%)' % (np.round(pca.var_exp[1] * 100)), fontsize=19)

    plt.title(figure_title, fontsize=18)
    plt.tight_layout()
    plt.savefig(figure_name, dpi=600)
    plt.close()


def merge_sample_and_metadata_annots(sample_annotation_df: pd.DataFrame,
                                     meta_data: pd.DataFrame,
                                     remove_qc_failed: bool,
                                     keep_replicates: bool = False,
                                     keep_reference: bool = False) -> pd.DataFrame:
    if remove_qc_failed:
        if 'QC' in sample_annotation_df.columns:
            sample_annotation_df = sample_annotation_df[sample_annotation_df['QC'] == 'passed']
        else:
            sample_annotation_df = sample_annotation_df[sample_annotation_df['Failed'] != 'x']

    meta_data_samples = meta_data['Sample name'].tolist()
    if keep_replicates:  # MT: simplify if statements
        sample_annotation_df = get_replicate_groups(sample_annotation_df, meta_data_samples)

    if not keep_replicates and keep_reference:
        sample_annotation_df = sample_annotation_df.loc[
                               (sample_annotation_df['Sample name'].isin(meta_data_samples)) | (
                                   sample_annotation_df['Sample name'].str.startswith('Reporter')), :]

    if not keep_replicates and not keep_reference:
        sample_annotation_df = meta_data.loc[meta_data['Sample name'].isin(sample_annotation_df['Sample name'].tolist()), :]

    # merge the two annotations to keep replicates+reference if any
    merged_annot = pd.merge(left=sample_annotation_df, right=meta_data, how='left', left_on='Sample name',
                            right_on='Sample name', suffixes=('', '_drop'))

    # drop columns present in both dataframes
    merged_annot = merged_annot.drop([col for col in merged_annot.columns if '_drop' in col], axis=1)

    # if null in merged annot retrieve batch from sample_annot
    for i, value in enumerate(merged_annot['Batch_No']):
        if math.isnan(value):
            batch = sample_annotation_df.loc[
                sample_annotation_df['Sample name'] == merged_annot.loc[i, 'Sample name'], 'Batch Name']
            merged_annot.loc[merged_annot.index.tolist()[i], 'Batch_No'] = batch.values[0]
    merged_annot['Batch_No'] = merged_annot['Batch_No'].astype(int)

    merged_annot['Tumor cell content'] = merged_annot['Tumor cell content'].apply(group_tumor_cell_content)
    merged_annot['Is reference channel'] = merged_annot['Sample name'].str.startswith('Reporter')
    merged_annot = merged_annot.replace(np.nan, 'nan')
    return merged_annot


def get_replicate_groups(sample_annotation_df, meta_data_samples):
    replicates = sample_annotation_df[
        ~(sample_annotation_df['Sample name'].isin(meta_data_samples)) &
        ~(sample_annotation_df['Sample name'].str.startswith('Reporter'))]

    # split off the replicate identifier from the sample name, e.g. ABCD-R2 => ABCD
    replicates['Sample name'] = replicates['Sample name'].apply(lambda x: '-'.join(x.split('-')[:-1]))
    sample_annotation_df['Replicate group'] = np.nan

    for i, sample in enumerate(replicates['Sample name'].unique()):  # with one there is the one without and one with -R2
        group_indices = sample_annotation_df[sample_annotation_df['Sample name'].str.contains(sample)].index
        sample_annotation_df.loc[group_indices.values, 'Replicate group'] = i + 1

    return sample_annotation_df


def group_tumor_cell_content(value) -> str:
    try:
        new_value = int(float(value))
    except ValueError:
        if pd.isnull(value):
            new_value = 'missing'
        elif '-' in value:
            the_range = value.split('-')
            range_average = sum([float(j) for j in the_range]) / len(the_range)
            new_value = round(range_average)
        elif any(string in str(value) for string in
                 ['missing', 'n.d.', '?', 'nd', 'ND (benigne)']):  # MT: check for this before trying to parse as float
            new_value = 'missing'
        else:
            logger.info('ERROR, what is this value: ', value)

    if new_value != 'missing':
        if new_value < 10:
            new_value = '0-10%'
        elif new_value < 20:
            new_value = '10-20%'
        elif new_value < 30:
            new_value = '20-30%'
        elif new_value < 40:
            new_value = '30-40%'
        elif new_value < 60:
            new_value = '40-60%'
        elif new_value < 80:
            new_value = '60-80%'
        else:
            new_value = '80-100%'
    return new_value


def create_ref_sample_annot(results_folder, sample_annot_df):
    df = pd.read_csv(os.path.join(results_folder, f'annot_fp_with_ref.csv'), index_col='Gene names')
    df = utils.keep_only_sample_columns(df)
    for sample in df:
        if sample not in sample_annot_df['Sample name'].tolist():
            # retrieve tmt channel, batch
            tmt_channel = re.search('corrected \d{1,2}', sample).group().split(' ')[1]
            batch = re.search('Batch\d{1,2}', sample).group()
            batch = re.findall(r'\d+', batch)[0]
            cohort = re.search('\d [A-Z,a-z]+_', sample).group().split(' ')[1][:-1]
            new_sample = {'Sample name': sample, 'Cohort': cohort, 'Batch Name': batch, 'TMT Channel': tmt_channel,
                          'QC': 'passed'}
            sample_annot_df = sample_annot_df.append(new_sample, ignore_index=True)
    return sample_annot_df


def whitespace_remover(df):
    for col in df.columns:
        if df[col].dtype == 'object':
            # applying strip function on column
            df[col] = df[col].astype(str).map(str.strip)
        else:
            pass
    return df


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of `x` and `y`

    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this two-dimensionl dataset
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = matplotlib.patches.Ellipse((0, 0),
                                         width=ell_radius_x * 2,
                                         height=ell_radius_y * 2,
                                         facecolor=facecolor,
                                         **kwargs)

    # Calculating the stdandard deviation of x from the squareroot of the variance and multiplying
    # with the given number of standard deviations
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = matplotlib.transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", dest="results", default=None, metavar="DIR", required=True, help="Path to results to check")

    argv = sys.argv[1:]
    args = parser.parse_args(argv)
    results_folder = Path(args.results)

    os.makedirs(os.path.join(results_folder, 'Results_investigation_plots'), mode=0o700, exist_ok=True)

    config_file = os.path.join(results_folder, 'configs.json')
    configs = config.load(config_file)
    # plot_types = ['Basket', 'Kinase', 'Intensity', 'Intensity_subset']
    # plot_types = ['Basket', 'Kinase']
    # plot_types = ['Kinase']
    plot_types = ['Intensity']
    # plot_types = ['Intensity_subset']
    # data_types = configs['data_types']
    data_types = ['pp']
    # create_plots(results_folder, plot_types, configs['sample_annotation'], configs['metadata_annotation'], configs["data_types"],
    #              min_sample_occurrence_ratio=0.5, umap=False, ref=False, replicate=True)

    # percentages = [0.5, 0.7, 1]
    percentages = [0.5]
    selected_proteins = []
    for min_sample_occurrence_ratio in percentages:
        for data_type in data_types:
            if data_type == 'pp':
                logger.info('data type: ', data_type)
                all_principal_dfs, all_principal_variances, meta_data_dicts = do_pca(selected_proteins, results_folder, plot_types,
                                                                                     configs['sample_annotation'],
                                                                                     configs['metadata_annotation'], data_type,
                                                                                     min_sample_occurrence_ratio,
                                                                                     umap=False, include_reference_channels=False,
                                                                                     include_replicates=True, plot=True)
