from __future__ import print_function
import os
import re

import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('pdf')

import warnings
import matplotlib.pyplot as plt
from job_pool import JobPool

from .. import mtb_scoring


def create_and_save_swarmplots_parallel(plot_function, plot_folder, plot_function_add_dot, df, num_threads, channels, *args):
    os.makedirs(plot_folder, mode=0o700, exist_ok=True)
    ax = create_boxplot(df)

    processingPool = JobPool(processes=num_threads)
    for channel in channels:
        processingPool.applyAsync(plot_function, [plot_folder, plot_function_add_dot, df, ax, channel, *args])

    processingPool.checkPool(printProgressEvery=num_threads)


def create_and_save_basket_boxplot(plot_folder, plot_function_add_dot, df, ax, channel, data_type, number_table, suffix):
    ax = plot_function_add_dot(df, ax, channel)
    title_suffix = ""
    if number_table is not None:
        temp_number = number_table[number_table['Sample'].str.startswith(channel)]
        protein_number = temp_number[temp_number['variable'] == f'Number proteins {data_type.lower()} rtk'][
            'value'].values
        title_suffix = "    " + 'n=' + str(protein_number[0])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel("Basket score", fontsize=16)
    ax.set(xlabel=None)
    title = channel + ', ' + data_type + title_suffix
    plt.title(title, fontsize=18)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)
    plt.tight_layout()
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning, message=r".*of the points cannot be placed.*")
        plt.savefig(os.path.join(plot_folder, f'{channel}_{data_type}_score_{suffix}.pdf'))
    plt.close()


def create_boxplot(df):
    plt.figure(figsize=(12, 8))
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=UserWarning, message=r".*of the points cannot be placed.*")
        ax = sns.boxplot(x="variable", y="value", data=df, color='lightgrey')
        ax = sns.stripplot(x="variable", y="value", data=df,
                           color=".25", size=2)
        plt.ylim([-150, 150])
    return ax


def add_patient_dot(df, ax, channel):
    temp_patient = df[df['Sample'].str.startswith(channel)]
    ax = sns.stripplot(x="variable", y="value", data=temp_patient,
                           color="red", size=2, ax=ax)
    return ax


def add_different_color_dots(df, ax, channel):
    all_colors = ['deepskyblue', 'darkorange', 'limegreen', 'mediumorchid']

    different_sample_types = df.loc[df['Sample'].str.startswith(channel), 'Sample'].unique()
    colors = all_colors[0:len(different_sample_types)]
    for i, var in enumerate(different_sample_types):
        temp_patient = df[df['Sample'].str.startswith(var)]
        ax = sns.stripplot(x="variable", y="value", data=temp_patient,
                               color=colors[i], size=2, ax=ax)
    return ax


def add_batch_color_dots(df, ax, channel):
    all_colors = ['lightgrey', 'dodgerblue', 'darkorange', 'limegreen', 'mediumorchid', 'silver', 'sienna', 'darkturquoise', 'darkkhaki',
                  'blueviolet', 'seagreen', 'orange', 'tan', 'olive', 'lightpink',
                  'mediumvioletred', 'khaki', 'powderblue', 'lightsalmon', 'olivedrab', 'firebrick',
                  'lawngreen', 'steelblue', 'indigo', 'linen', 'springgreen', 'gold', 'darkred', 'lightgreen', 'pink',
                  'yellowgreen', 'red', 'blue', 'lightsteelblue', 'palegoldenrod', 'coral', 'yellow', 'magenta',
                  'darkcyan', 'maroon', 'peru', 'palevioletred', 'navy', 'orange', 'bisque', 'orangered', 'rosybrown', 'gold',
                  'slateblue', 'palegreen', 'darkseagreen', 'deepskyblue', 'crimson', 'deeppink', 'mediumvioletred']
    different_sample_types = df.loc[df['Sample'].str.startswith(channel), 'Batch'].unique()
    colors = all_colors[0:len(different_sample_types)]

    for i, var in enumerate(different_sample_types):
        temp_patient = df[df['Batch'] == var]
        ax = sns.stripplot(x="variable", y="value", data=temp_patient,
                               color=colors[i], size=2, ax=ax)
    return ax


def get_ref_outliers(df):
    all_outliers = []
    for basket in df.groupby('variable'):
        basket = basket[1]
        print(basket.head(1))
        ref = basket[basket['Sample'].str.startswith('Reporter intensity')]
        print(ref.value.median())
        print(ref.value.std())
        outliers = ref.loc[(ref['value'] > ref.value.median()+abs(2*ref.value.std())) | (ref['value'] < ref.value.median()-abs(2*ref.value.std())), :]
        all_outliers.append(outliers)
    all_outliers = pd.concat(all_outliers, axis=0)
    return all_outliers


if __name__ == "__main__":
    zug_office = True
    folder = '2022.10.26_CJ_batch57'
    folder = '2022.11.18_CJ_batch59_test'
    results_folder = f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}'
    # results_folder = f'/home/cjensen/Desktop/test_data_plot/'
    basket_scores = pd.read_csv(f'{results_folder}/basket_scores_ref.tsv', sep='\t', index_col='Sample')

    clinical_baskets, immune_rtk_baskets, rtk_baskets = mtb_scoring.get_baskets_order(basket_scores.columns)
    rtk_number_table = mtb_scoring.create_num_rtk_annotation_df(basket_scores, data_types='fp')
    rtk_number_table = rtk_number_table[rtk_number_table['Sample'].str.startswith('Reporter intensity')]
    basket_types = [('baskets_ref', 'Scores', clinical_baskets, None),
                    ('immune_rtk', 'Scores', immune_rtk_baskets, None),
                    ('rtk_baskets_ref', 'RTK Scores', rtk_baskets, rtk_number_table)]

    # bool if we want to subset some batches  - this is more for the swarm plot...
    subset_batches = True

    for data_type in ['FP', 'PP']:
        if data_type == 'PP':
            if 'Immunotherapy' in clinical_baskets:
                clinical_baskets.remove('Immunotherapy')
                clinical_baskets.remove('DRD')

        for basket_type, sheet_suffix, baskets, number_table in basket_types:
            if len(baskets) == 0:
                print(f"Skipping {basket_type} {data_type} box plots as no baskets were given")
                continue

            print(f"Plotting {basket_type} {data_type} box plots")
            df = mtb_scoring.create_boxplot_df(basket_scores, f'{data_type} - {sheet_suffix}', baskets)

            # get tables for outliers
            outliers = get_ref_outliers(df)
            # save outliers
            outliers.to_csv(os.path.join(results_folder, f'{data_type}_{basket_type}_channel_outliers.csv'))
            channels = ['Reporter intensity corrected 9', 'Reporter intensity corrected 10', 'Reporter intensity corrected 11']
            df['Batch'] = df['Sample']
            df['Sample'] = df['Batch']

            for channel in channels:
                for i, ele in enumerate(df['Sample']):
                    if channel in ele:
                        test_reg = re.search(r'[a-z,A-Z]+_Batch\d{1,2}', ele)
                        df.loc[i, 'Sample'] = channel
                        df.loc[i, 'Batch'] = test_reg.group(0)
            all_samples = df['Sample']
            num_threads = 1

            # TODO make the swarmplot stay the same and just the coloring be different??

            plot_folder = os.path.join(results_folder, basket_type+'_all_channels_together')
            create_and_save_swarmplots_parallel(create_and_save_basket_boxplot, plot_folder, add_patient_dot, df, num_threads,
                                                ['Reporter intensity corrected'], data_type, number_table, basket_type)

            plot_folder = os.path.join(results_folder, basket_type+'_each_channel_alone')
            create_and_save_swarmplots_parallel(create_and_save_basket_boxplot, plot_folder, add_patient_dot, df, num_threads,
                                                channels, data_type, number_table, basket_type)

            plot_folder = os.path.join(results_folder, basket_type+'_all_channels_together_each_their_color')
            create_and_save_swarmplots_parallel(create_and_save_basket_boxplot, plot_folder, add_different_color_dots, df,
                                                num_threads, ['Reporter intensity corrected'], data_type, number_table, basket_type)

            plot_folder = os.path.join(results_folder, basket_type+'_all_channels_together_color_per_batch')
            create_and_save_swarmplots_parallel(create_and_save_basket_boxplot, plot_folder, add_batch_color_dots, df,
                                                num_threads, ['Reporter intensity corrected'], data_type, number_table, basket_type)
