import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


#folder_1 = '2022.08.12_CJ_debug_mode_minimal_no_simsi'
#folder_2 = '2022.08.12_CJ_debug_mode_minimal_simsi'
folder_1 = '2023.09.29_CJ_minimal_test_plotting'

# choose 3-5 batches (first ones)
data_type = 'fp'

# Both datasets
# before_median_centering = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_fp2_before_ms1_correction.csv', index_col='Gene names', nrows=100)
# before_median_centering = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preproces'
#                                       f'sed_{data_type}_complete_raw.csv', index_col='Gene names')
# print(before_median_centering.head())

# after_median_centering = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_{data_type}_after_1st_median.csv', index_col='Gene names')
# TODO: add before ms1 median and make them go into one boxplot
after_ms1_median_centering = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_{data_type}_after_ms1_centering.csv', index_col='Gene names')
# after_ms1_imputation = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_{data_type}_after_ms1_imputation.csv', index_col='Gene names', nrows=100)
# after_ms1_redistribution = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_{data_type}_after_ms1_correction.csv', index_col='Gene names', nrows=100)
# before_removing_qc_failed = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_{data_type}_before_filter_sample.csv', index_col='Gene names', nrows=100)

# Phospho only
# after_imputation = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_pp_after_imputation.csv', index_col='Gene names', nrows=100)
# after_sum_peptide = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/debug_preprocessed_pp_after_aggregation.csv',
#                                 index_col='Gene names', nrows=100)

# dfs = [before_median_centering, after_median_centering]
dfs = [after_ms1_median_centering]
# df_labels = ['raw', 'after_median', 'after_ms1_median', 'after_ms1_imp', 'after_ms1_redist', 'before_qc_filter']
# df_labels = ['raw', 'after_median']
df_labels = ['after_ms1_median']

for i, df in enumerate(dfs):

    # full proteome

    if 'Batch' in df.columns:
        df = df.reset_index()
        df = df.set_index(['Gene names', 'Batch'])
    else:
        df = df.reset_index()
        df = df.set_index(['Gene names', 'Experiment'])

    if df_labels[i] not in ['after_ms1_median', 'after_ms1_imp']:
        if df.columns.str.contains('Reporter intensity corrected').any():
            df = df.filter(regex='Reporter intensity corrected \d{1,2}')
        else:
            df = df.filter(regex=r'^\S+-.+-\S+')
    else:
        # subset to MS1 instead
        df = df.filter(regex='^Intensity')

    print(df.head())
    # if needed
    df = np.log10(df)
    # subset to batch 1-3  - todo: change to random ones
    df = df.reset_index()
    df = df.set_index('Gene names')
    batches = df['Batch'].unique()[0:3]
    df = df[df['Batch'].isin(batches)]
    print(df.shape)
    # batch 1, 27, 54 and else as above

    # print(df.head())
    if df_labels[i] != 'after_ms1_median':
        df = df.groupby('Batch')
    # all_batches, all_batches_names = [], []
    # for name, batch in df:
    #     batch = batch.reset_index()
    #     batch = batch.set_index('Gene names')
    #     batch = batch.drop(columns='Batch')
    #     all_batches.append(batch)
    #     all_batches_names.append(name)
    #
    # new_df = pd.merge(all_batches[0], all_batches[1], left_index=True, right_index=True)
    # print(new_df.shape)
    # print(new_df.head())
    #
    # exit()
    # new_df = pd.merge(new_df, all_batches[2], left_index=True, right_index=True)
    # print(new_df.shape)

    # for batch in df:
    for name, batch in df:
        fig = plt.figure(figsize=(12, 10))
        sns.boxplot(batch)
        plt.xticks(rotation=45)
        plt.ylim(-2, 12)
        plt.tight_layout()
        plt.savefig(f'/home/cjensen/kusterlab_home/kusterlab/users_files/Cecilia_Jensen/Preprocessing/boxplot_{data_type}_{df_labels[i]}_{name}.png')
        # plt.savefig(f'/media/kusterlab/users_files/Cecilia_Jensen/Preprocessing/boxplot_{data_type}_{df_labels[i]}_{name}.png')

# per data prepare
# sns.boxplot(data=df, )



# plot





exit()

# df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/preprocessed_fp2_before_ms1_correction.csv', index_col='Gene names')
# df_filtered = df.filter(regex='Reporter intensity corrected \d{1,2}')
# df_filtered['Batch'] = df['Batch']
# df_filtered = df_filtered.groupby('Batch')
# df = df.filter(regex=r'^\S+-.+-\S+')

# df2 = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_2}/preprocessed_fp2_after_ms1_correction.csv', index_col='Gene names')
# # df2 = df2.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+')
# df2_filtered = df2.filter(regex='Reporter intensity corrected \d{1,2}')
# df2_filtered['Batch'] = df2['Batch']
# df2_filtered = df2_filtered.groupby('Batch')


# boxplot
# for batch, df in df_filtered:
#     print(batch)
#     # take each batch separately and subtract the mean ref channel
#     # intensity_batch = df.iloc[:, 11 * i:(11 * (i + 1))]
#     # intensity_batch = intensity_batch.iloc[:, [0, 1, 2, 3, 4, 5, 6, 7]]
#     intensity_batch = df.replace(0, np.nan)
#     intensity_batch = intensity_batch.filter(regex='Reporter intensity corrected \d{1,2}')
#     # # if corrected is False:
#     intensity_batch = np.log10(intensity_batch)
    # print(intensity_batch)

    # Add title with batch # and set limits of axis to a constant

    # # create a list of length values
    # batch1_labels = ['sample 1' for i in range(len(intensity_batch.index))]
    # batch2_labels = ['sample 2' for i in range(len(intensity_batch.index))]
    # batch3_labels = ['sample 3' for i in range(len(intensity_batch.index))]
    # batch4_labels = ['sample 4' for i in range(len(intensity_batch.index))]
    # batch5_labels = ['sample 5' for i in range(len(intensity_batch.index))]
    # batch6_labels = ['sample 6' for i in range(len(intensity_batch.index))]
    # batch7_labels = ['sample 7' for i in range(len(intensity_batch.index))]
    # batch8_labels = ['sample 8' for i in range(len(intensity_batch.index))]
    # all_batch_labels = [batch1_labels, batch2_labels, batch3_labels, batch4_labels, batch5_labels, batch6_labels,
    #                     batch7_labels, batch8_labels]
    # all_batch_labels = [element for sublist in all_batch_labels for element in sublist]
    # all_columns = [intensity_batch.iloc[:, 0], intensity_batch.iloc[:, 1], intensity_batch.iloc[:, 2], intensity_batch.iloc[:, 3],
    #                intensity_batch.iloc[:, 4], intensity_batch.iloc[:, 5], intensity_batch.iloc[:, 6], intensity_batch.iloc[:, 7]]
    # all_columns = pd.concat(all_columns)
    # # make df and add sample column
    # all_columns = pd.DataFrame({'Intensities': all_columns,
    #                             'Samples': all_batch_labels})
    # print(all_columns)
    #
    # f, ax = plt.subplots(figsize=(12, 10))
    # sns.boxplot(x='Samples', y='Intensities', data=all_columns, color='grey')
    # ax.set_title(f'{batch}', fontsize=22)
    # ax.set_xlabel('Samples', fontsize=22)
    # ax.set_ylabel('Intensities', fontsize=22)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=22)
    # # ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize=22)
    # ax.tick_params(labelsize=22)
    # plt.ylim(0, 12)

    # # save plots
    # # if corrected:
    # #     plt.savefig(
    # #         '/home/cjensen/kusterlab/internal_projects/active/TOPAS/WP31'
    # #         '/Searches/Patients/Analysis_CJ/boxplots/PP_Batch'
    # #         '%s_normalized'
    # #         '.png' % str(
    # #             i + 1))
    # #     plt.close()
    # # else:
    # plt.savefig(
    #     f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_2}/before_ms1_correction_{batch}.jpg', bbox_inches='tight', dpi=400)
    # plt.close()


