import os
import re
import json
import math
from typing import List, Dict, Tuple, Any, Union
from itertools import product

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from ppca import PPCA

from umap import UMAP

from .. import check_results_quality

# def read_data():


if __name__ == "__main__":
    # Parameters
    zug_office = True
    # folder = '2022.10.25_CJ_minimal_test_ref'
    folder = '2022.10.26_CJ_batch57'
    # results_folder = f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}'
    results_folder = f'/home/cjensen/Desktop/test_data_plot/'
    config_file = os.path.join(results_folder, 'configs.json')

    invest_type = 'basket_score'  # one of ['all', 'baskets', 'basket_score']
    pct_threshold_identificated = 0.8  # if invest_type all is chosen this here is used for filtering before imputation
    targets_name = 'Entities'  # getting targets and title of legend

    # Make different umaps.. loop through parameters
    # num_components = [10, 30, 50]
    num_components = 30
    # num_neighbors = [5, 10, 15, 30]
    num_neighbors = 5
    # metric = ['euclidean', ]
    metric = 'euclidean'
    # epochs = [100, 1000, 10000]
    epochs = 1000

    # TODO new label of drug recommendation?

    with open(config_file, "r") as inp:
        configs = json.load(inp)
    if zug_office:
        configs = {}
        configs['sample_annotation'] = f'{results_folder}/patient_annotation_mixed_cohort_221026.csv'
        configs['metadata_annotation'] = f'{results_folder}/MasterInform_Metadata_221014_CJ.xlsx'

    # Read in annotation files
    sample_annotation_df = pd.read_csv(configs['sample_annotation'])
    metadata_df = pd.read_excel(configs['metadata_annotation'])
    metadata_df['Sample name'] = metadata_df['Sample name'].str.strip()

    # metadata_columns = ['Replicates', 'Batch_No', 'Sarcoma Subtype', 'Program', 'Tumor cell content']
    # metadata_dicts = {}
    # for metadata in metadata_columns:
    #     metadata_dicts[metadata] = check_results_quality.get_sample_dict(sample_annotation_df, metadata_df, remove_qc_failed=True,
    #                                                                      data_type=metadata)

    if invest_type in ['all', 'baskets']:
        # Read in intensities
        # intensity_df = pd.read_csv(f'{results_folder}/preprocessed_fp.csv', index_col='Gene names')
        # df = intensity_df.loc[:, ~intensity_df.columns.str.startswith('Identification')]
        intensity_df = pd.read_csv(f'{results_folder}/annot_fp.csv', index_col='Gene names')
        annot = intensity_df[['basket', 'rtk']]
        df = intensity_df.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+')

    if invest_type == 'baskets':
        # Subset to basket annotation
        main_bool = annot['basket'].isna()
        rtk_bool = annot['rtk'].isna()
        df_main = df.loc[main_bool, :]
        df_rtk = df.loc[rtk_bool, :]

    if invest_type == 'basket_score':
        basket_scores = pd.read_csv(f'{results_folder}/basket_scores.tsv', sep='\t', index_col='Sample')

        # separate fp and pp
        basket_scores_fp = basket_scores.iloc[:, basket_scores.columns.str.startswith('FP')]
        basket_scores_pp = basket_scores.iloc[:, basket_scores.columns.str.startswith('PP')]
        # remove drug scores
        basket_scores_fp = basket_scores_fp.iloc[:, ~basket_scores_fp.columns.str.contains('TOPAS')]
        basket_scores_pp = basket_scores_pp.iloc[:, ~basket_scores_pp.columns.str.contains('TOPAS')]
        # remove RTK baskets
        basket_scores_fp_main = basket_scores_fp.iloc[:, ~basket_scores_fp.columns.str.contains('RTK')]
        basket_scores_pp_main = basket_scores_pp.iloc[:, ~basket_scores_pp.columns.str.contains('RTK')]

    df = basket_scores_fp_main  # df_main
    df = df.transpose()

    if targets_name == 'Batches':
        targets = []
        for sample in df.columns:
            targets.append(sample_annotation_df.loc[sample_annotation_df['Sample name'] == sample, 'Batch Name'].values[0])

    if targets_name == 'Entities':
        targets = []
        for sample in df.columns:
            # remove -R if any
            patient_match = re.search('[A-Z]\d+-[A-Z|\d]+-[A-Z|\d]+-R\d', sample)
            if patient_match:
                sample = '-'.join(sample.split('-')[:-1])

            if sample in metadata_df['Sample name'].tolist():
                targets.append(metadata_df.loc[metadata_df['Sample name'] == sample, 'Sarcoma Subtype'].values[0])

    # Transpose and filter for threshold
    df = df.transpose()
    df = df.loc[:, df.count(axis=0) >= df.shape[0] * 0.8]

    # Calculate imputed data with PPCA
    ppca = PPCA()
    x = StandardScaler().fit_transform(df)
    ppca.fit(data=x, d=2, verbose=False)

    df1 = pd.DataFrame(x)
    imputed_data = pd.DataFrame(ppca.data)
    imputed_data.index = df.index
    imputed_data.columns = df.columns



    # Calculate UMAP
    reducer = UMAP(n_components=components, n_epochs=epochs, low_memory=True, n_neighbors=neighbors, metric=metric)
    # reducer = UMAP(n_components=30, n_epochs=1000, low_memory=True, n_neighbors=5, metric='euclidean')
    X_trans = reducer.fit_transform(imputed_data)

    # Prepare dataframe
    umap_df = pd.DataFrame(data=X_trans)
    umap_df = pd.concat([umap_df.iloc[:, 2], umap_df.iloc[:, 3], pd.Series(df.index), pd.Series(targets)], axis=1)
    umap_df.columns = ['Principal component 1', 'Principal component 2', 'Sample', 'Targets']

    # TODO prepare different labels

    # PLOT
    all_colors = ['lightgrey', 'dodgerblue', 'darkorange', 'limegreen', 'mediumorchid', 'mediumvioletred', 'sienna',
                  'darkturquoise', 'darkkhaki', 'blueviolet', 'seagreen', 'orange', 'slategrey', 'tan', 'olive', 'lightpink',
                  'darkslategray', 'mediumvioletred', 'khaki', 'powderblue', 'lightsalmon', 'olivedrab', 'firebrick',
                  'lawngreen', 'steelblue', 'indigo', 'linen', 'springgreen', 'gold', 'darkred', 'lightgreen', 'pink',
                  'yellowgreen', 'red', 'blue', 'lightsteelblue', 'palegoldenrod', 'coral', 'yellow', 'magenta', 'darkslategrey',
                  'darkcyan', 'maroon', 'peru', 'palevioletred', 'navy', 'orange', 'bisque', 'orangered', 'rosybrown', 'gold',
                  'slateblue', 'palegreen', 'darkseagreen', 'deepskyblue', 'crimson', 'deeppink']

    colors = all_colors[0:len(umap_df['Targets'].unique())]
    pca_combi = zip(umap_df['Targets'].unique(), colors)
    fig, ax = plt.subplots(1, 1, figsize=(8.27, 8.27))
    for target, color in pca_combi:
        temp_df = umap_df.loc[umap_df['Targets'] == target, :]
        ax.scatter(temp_df.loc[:, 'Principal component 1'],
                   temp_df.loc[:, 'Principal component 2'], s=25, c=color)

    legend_size = 6
    if len(umap_df['Targets'].unique()) > 30:
        legend_size = 5
    ax.legend(umap_df['Targets'].unique(), prop={'size': legend_size}, loc='center right', bbox_to_anchor=(1.15, 0.5),
              title=f'{targets_name}')
    plt.title(f'UMAP of Sarcoma cohort - {invest_type}')
    plt.tight_layout()
    plt.savefig(
        os.path.join(results_folder, 'Results_investigation_plots', f'Umap_{invest_type}_{num_neighbors}_{num_components}_{metric}2.pdf'))

    # basket score data entity
