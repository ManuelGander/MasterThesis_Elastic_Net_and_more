import os
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

from .. import check_results_quality


def calculate_pcas(df: pd.DataFrame,
                   pct_threshold: List[float] = [0.5],
                   pca: bool = False) -> Tuple[List, List]:
    """Calculate PCA or PPCA of given df with given threshold"""
    all_x, pcas = [], []
    for i, pct in enumerate(pct_threshold):
        # define dataframe with threshold and transform it
        x = df.loc[:, df.count(axis=0) >= df.shape[0] * pct]
        x = StandardScaler().fit_transform(x)
        all_x.append(x)
        # if do_pca is true do first an actual pca
        if pca and i == 0:
            pca = PCA(n_components=2)
            pca.principal_components = pca.fit_transform(x)
            pcas.append(pca)
            ppca = PPCA()
            ppca.fit(data=x, d=2, verbose=False)
            ppca.var_exp[1] = ppca.var_exp[1] - ppca.var_exp[0]
        else:
            ppca = PPCA()
            ppca.fit(data=x, d=2, verbose=False)
            # change cumulated explained variance to just variance per PC
            ppca.var_exp[1] = ppca.var_exp[1] - ppca.var_exp[0]
        pcas.append(ppca)
    if len(pcas) == 1:
        pcas = pcas[0]
    return pcas, all_x


if __name__ == "__main__":
    # folder = '2022.10.25_CJ_minimal_test_ref'
    zug_office = False
    folder = '2022.11.18_CJ_batch59_test'
    results_folder = f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}'
    annot_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches'
    meta_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation'
    if zug_office:
        folder = '2022.10.26_CJ_batch57'
        configs = {}
        configs['sample_annotation'] = f'{results_folder}/patient_annotation_mixed_cohort_221026.csv'

    # results_folder = f'/home/cjensen/Desktop/test_data_plot/'
    config_file = os.path.join(results_folder, 'configs.json')
    with open(config_file, "r") as inp:
        configs = json.load(inp)
    configs['sample_annotation'] = f'{annot_folder}/patient_annotation_mixed_cohort_221031.csv'
    configs['metadata_annotation'] = f'{meta_folder}/MasterInform_Metadata_221014_CJ.xlsx'

    # Read in annotation files
    sample_annotation_df = pd.read_csv(configs['sample_annotation'])
    metadata_df = pd.read_excel(configs['metadata_annotation'])
    metadata_df['Sample name'] = metadata_df['Sample name'].str.strip()

    # Read in basket scores # TODO: make it also without ref and then loop over
    basket_scores = pd.read_csv(f'{results_folder}/basket_scores_ref.tsv', sep='\t', index_col='Sample')
    # basket_scores = pd.read_csv(f'{results_folder}/basket_scores.tsv', sep='\t', index_col='Sample')
    # separate fp and pp
    basket_scores_fp = basket_scores.iloc[:, basket_scores.columns.str.startswith('FP')]
    basket_scores_pp = basket_scores.iloc[:, basket_scores.columns.str.startswith('PP')]
    # remove drug scores
    basket_scores_fp = basket_scores_fp.iloc[:, ~basket_scores_fp.columns.str.contains('TOPAS')]
    basket_scores_pp = basket_scores_pp.iloc[:, ~basket_scores_pp.columns.str.contains('TOPAS')]
    # remove RTK baskets
    basket_scores_fp_main = basket_scores_fp.iloc[:, ~basket_scores_fp.columns.str.contains('RTK')]
    basket_scores_pp_main = basket_scores_pp.iloc[:, ~basket_scores_pp.columns.str.contains('RTK')]

    # Create data dictionary
    data_dict = {'preprocessed_fp_ref': basket_scores_fp_main, 'preprocessed_pp_ref': basket_scores_pp_main}

    # metadata_columns = ['Replicates', 'Batch_No', 'Sarcoma Subtype', 'Program', 'Tumor cell content']
    # metadata_dicts = {}
    # for metadata in metadata_columns:
    #     metadata_dicts[metadata] = check_results_quality.get_sample_dict(sample_annotation_df, metadata_df, remove_qc_failed=True,
    #                                                                      data_type=metadata)

    # TODO have a script with constants such as plot colors?
    # # TODO: color based on batch, replicates,
    all_colors = ['lightgrey', 'dodgerblue', 'darkorange', 'limegreen', 'mediumorchid', 'silver', 'sienna', 'darkturquoise', 'darkkhaki',
                  'blueviolet', 'seagreen', 'orange', 'slategrey', 'tan', 'olive', 'lightpink', 'darkslategray',
                  'mediumvioletred', 'khaki', 'powderblue', 'lightsalmon', 'olivedrab', 'firebrick',
                  'lawngreen', 'steelblue', 'indigo', 'linen', 'springgreen', 'gold', 'darkred', 'lightgreen', 'pink',
                  'yellowgreen', 'red', 'blue', 'lightsteelblue', 'palegoldenrod', 'coral', 'yellow', 'magenta', 'darkslategrey',
                  'darkcyan', 'maroon', 'peru', 'palevioletred', 'navy', 'orange', 'bisque', 'orangered', 'rosybrown', 'gold',
                  'slateblue', 'palegreen', 'darkseagreen', 'deepskyblue', 'crimson', 'deeppink', 'mediumvioletred']
    all_symbols = ['o', 'v', '^', 's', '+', 'd', 'x', 'p', 'P', '<', '>']

    principal_dfs = []
    data_type = ['fp', 'pp']
    for i, df in enumerate(data_dict):
        df = data_dict[df]
        # print(df)

        # TODO: make it a function to retrieve either replicates or batch annot as targets
        # how many replicates
        is_replicate = pd.Series(df.index.str.split('-').str[-1])
        is_replicate = is_replicate.apply(lambda x: x if 'R' in x and 'Reporter' not in x else '')
        sample_no_replicate = pd.Series(df.index.str.replace(r'-R[0-9]', ''))
        is_replicate = pd.concat([is_replicate, sample_no_replicate, pd.Series(df.index)], axis=1)
        is_replicate.columns = ['Replicate', 'Sample', 'Samples']
        # sample_targets = is_replicate.loc[is_replicate['Replicate'] != '', 'Sample'].unique()
        is_replicate['Replicate group'] = np.nan

        j = 1
        for sample in is_replicate['Sample'].unique():  # with one there is the one without and one with -R2
            group_indixes = is_replicate.loc[is_replicate['Samples'].str.contains(sample), :].index
            is_replicate.loc[group_indixes.values, 'Replicate group'] = j
            j += 1
        # Keep only replicate group numbers
        count_occurrences = is_replicate['Replicate group'].value_counts()
        count_occurrences = pd.Series(count_occurrences[count_occurrences == 1].index)
        is_replicate.loc[is_replicate['Replicate group'].isin(count_occurrences), 'Replicate group'] = ''

        targets = is_replicate['Replicate group']

        num_combinations = math.ceil(np.sqrt(len(targets.unique())))
        colors = all_colors[0:num_combinations]
        symbols = all_symbols[0:num_combinations]

        color_marker_tuples = []
        for shape, color in product(colors, symbols):
            color_marker_tuples.append((shape, color))
        color_marker_tuples = color_marker_tuples[0:len(targets.unique())]

        samples = pd.Series(df.index)
        samples = samples.apply(lambda x: x.replace('Reporter intensity corrected ', '') if 'Reporter' in x else '')
        batches = [sample.split('_')[1] if '_' in sample else '' for sample in samples]
        batches = pd.Series(batches)
        print(batches)
        color_dict = dict(zip(batches.unique(), all_colors[0:len(batches.unique())]))
        colors = [color_dict[batch] for batch in batches]

        samples = samples.str.get(0)+samples.str.get(1).fillna('')
        samples = samples.fillna('')

        targets = samples
        # print(targets.unique())


        pca_object, x = calculate_pcas(df, pct_threshold=[1], pca=True)

        # explained variance
        # TODO use this
        print(pca_object[1].var_exp)

        principal_df = pd.DataFrame(data=pca_object[1].transform(),
                                    columns=['Principal component 1', 'Principal component 2'])
        principal_df = pd.concat([principal_df, samples, targets], axis=1)
        principal_df.columns = ['Principal component 1', 'Principal component 2', 'Sample', 'Targets']
        principal_dfs.append(principal_df)

        fig, ax = plt.subplots(1, 1, figsize=(8.27, 8.27))
        # zip color+marker together with target
        pca_combi = zip(targets.unique(), color_marker_tuples)

        # for target, plot_tuple in pca_combi:
        #     temp_df = principal_df[principal_df['Targets'] == target]
        #     ax.scatter(temp_df.loc[:, 'Principal component 1'], temp_df.loc[:, 'Principal component 2'], s=30, c=plot_tuple[0], marker=plot_tuple[1])
        # legend_size = 6
        # if len(targets.unique()) > 30:
        #     legend_size = 4
        # ax.legend(targets.unique(), prop={'size': legend_size}, loc='center right', bbox_to_anchor=(1.15, 0.5))
        #
        # plt.tight_layout()
        # plt.savefig(os.path.join(results_folder, 'Results_investigation_plots', f'{data_type[i]}_Baskets_pca_ref_replicates.pdf'))
        # plt.savefig(os.path.join(results_folder, 'Results_investigation_plots', f'{data_type[i]}_Baskets_pca_replicates.pdf'))

    fig, axes = plt.subplots(2, 1, figsize=(8.27, 11.69))
    axes[0].scatter(principal_dfs[0].loc[:, 'Principal component 1'],
                    principal_dfs[0].loc[:, 'Principal component 2'], s=30, alpha=1, c=colors)
    axes[1].scatter(principal_dfs[1].loc[:, 'Principal component 1'],
                    principal_dfs[1].loc[:, 'Principal component 2'], s=30, alpha=1, c=colors)

    for i, txt in enumerate(principal_dfs[0]['Sample']):
        axes[0].annotate(txt, (principal_dfs[0].loc[i, 'Principal component 1'], principal_dfs[0].loc[i, 'Principal component 2']),
                         fontsize=8)
        axes[1].annotate(txt, (principal_dfs[1].loc[i, 'Principal component 1'], principal_dfs[1].loc[i, 'Principal component 2']),
                         fontsize=8)

    # TODO add legend
    # print(principal_dfs[0])
    # axes[0].legend()


    plt.savefig(os.path.join(results_folder, 'Results_investigation_plots', 'Baskets_pca_ref.pdf'))
