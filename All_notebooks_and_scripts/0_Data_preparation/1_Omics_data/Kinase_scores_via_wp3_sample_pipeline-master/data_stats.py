import os
import re
import sys
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

import bin.preprocess_tools as prep


def analyse_results(results_folder, sample_annotation, sample_for_report):
    """

    :param results_folder: location to data to analyse
    :param sample_annotation: not in use
    :param sample_for_report: not in use
    """

    try:
        # how about raw data??

        # fp_preprocessed_df = pd.read_csv('{0}/preprocessed_fp.csv'.format(
        #     results_folder))
        # fp_preprocessed_df.set_index('Protein IDs', inplace=True)
        # pp_preprocessed_df = pd.read_csv('{0}/preprocessed_pp.csv'.format(
        #     results_folder))
        # pp_preprocessed_df.set_index(['Modified sequence', 'Proteins', 'Gene names'],
        #                              inplace=True)
        #
        fp_clinically_processed = pd.read_csv('{0}/filtered_fp.csv'.format(
            results_folder))
        fp_clinically_processed.set_index('Protein IDs', inplace=True)
        fp_clinically_processed_subset = pd.read_csv('{0}/filtered_subset_fp.csv'
                                                     .format(results_folder))
        fp_clinically_processed_subset.set_index('Gene names', inplace=True)
        pp_clinically_processed = pd.read_csv('{0}/filtered_annot_pp.csv'.format(
            results_folder))
        pp_clinically_processed.set_index(
            ['Modified sequence', 'Proteins', 'Gene names'],
            inplace=True)
        pp_clinically_processed_subset = pd.read_csv('{0}/filtered_annot_subset_pp.csv'
                                                     .format(results_folder))
        pp_clinically_processed_subset.set_index(['Modified sequence', 'Proteins',
                                                  'Gene names'], inplace=True)

    except OSError as e:  # todo error doesnt happen when file is not there
        print(e)
        sys.exit(1)
    # print(fp_preprocessed_df.sum(axis=1))
    # print(fp_clinically_processed.sum(axis=1))
    # print(pp_preprocessed_df.iloc[:, 22:30])
    # temp_df = pp_preprocessed_df.iloc[:, 22:30]
    # # number of rows that are not missing
    # print('phosphopeptides in batch', temp_df.shape[0] - (temp_df.shape[0] -
    #                                                       temp_df.dropna().shape[0]))
    #
    # print('protein counts')
    #
    # # print(fp_preprocessed_df)
    # # temp_df = fp_preprocessed_df.iloc[:, 0:8]  # pp_preprocessed_df.iloc[:, 22:30]
    # # print(temp_df, temp_df.shape[0], temp_df.dropna().shape[0])
    # # number of rows that are not missing
    # # print('proteins in batch', temp_df.dropna().shape[0])
    # # just do it per column and ignore the ref channels
    # temp_df = fp_preprocessed_df.copy()
    # if len(temp_df.columns) % 11 == 0:
    #     n_batches = int(len(temp_df.columns) / 11)
    # else:
    #     print('Error number of columns is {}'.format(len(temp_df.columns)))
    #
    # ref_cols = [np.array([8, 9, 10]) + (11 * i) for i in range(n_batches)]
    # # flatten list of arrays to list
    # ref_cols = np.concatenate(ref_cols).ravel().tolist()
    # temp_df.drop(temp_df.columns[ref_cols], axis=1, inplace=True)
    # print(len(temp_df.columns))
    # start_index = [0]
    # i = 1
    # while i < 176/8:
    #     start_index.append(8*i)
    #     i += 1
    #
    # protein_counts = [0]
    # for i in start_index:
    #     print(temp_df.iloc[:, i:i+8])
    #     batch = temp_df.iloc[:, i:i+8].dropna()
    #     protein_counts.append(batch.shape[0])

    # print(protein_counts)
    # protein_counts = [7830, 6688, 6985, 4932, 6219, 5651, 5212, 6096, 5999, 5893, 6345,
    #                   4445, 7703, 7957, 7537, 6739, 6620, 8085, 8068, 7524, 7760, 7320]

    # plot
    #
    # peptide_counts = [22575, 20959, 21567, 13111, 17710, 18256, 12951, 16132, 16828,
    #                   13894, 11972, 21437, 20588, 20619, 20397, 17826, 23625, 25702,
    #                   26711, 20977, 19109, 21871]
    #
    # df = pd.DataFrame({
    #     'Factor': ['1', '2', '3', '4', '5',
    #                '6', '7', '8', '9', '10',
    #                '11', '16', '17', '18', '20',
    #                '21', '22', '23', '24', '25',
    #                '26', '27'],
    #     'Full proteome': protein_counts
    # })
    # fig, ax1 = plt.subplots(figsize=(11, 5))
    # tidy = df.melt(id_vars=['Factor']).rename(columns=str.title)
    # tidy.columns = ['Batch', 'Data', 'Number of proteins']
    # sns.barplot(x='Batch', y='Number of proteins', hue='Data', data=tidy, ax=ax1)
    # ax1.set_xlabel('Batches', fontsize=18)
    # ax1.set_ylabel('Proteins', fontsize=18)
    # ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=18)
    # plt.title('Full proteome', fontsize=22)
    # plt.legend([], [], frameon=False)
    # plt.tight_layout()
    # sns.despine(fig)
    # plt.savefig('/home/cjensen/Desktop/full_proteome_protein_count_per_batch.png')
    #
    # df = pd.DataFrame({
    #     'Factor': ['1', '2', '3', '4', '5',
    #                '6', '7', '8', '9', '10',
    #                '11', '16', '17', '18', '20',
    #                '21', '22', '23', '24', '25',
    #                '26', '27'],
    #     'Phospho proteome': peptide_counts
    # })
    # fig, ax1 = plt.subplots(figsize=(11, 5))
    # tidy = df.melt(id_vars=['Factor']).rename(columns=str.title)
    # tidy.columns = ['Batch', 'Data', 'Number of peptides']
    # sns.barplot(x='Batch', y='Number of peptides', hue='Data', data=tidy, ax=ax1)
    # ax1.set_xlabel('Batches', fontsize=18)
    # ax1.set_ylabel('Phospho-peptides', fontsize=18)
    # ax1.set_xticklabels(ax1.get_xticklabels(), fontsize=18)
    # plt.title('Phospho proteome', fontsize=22)
    # plt.legend([], [], frameon=False)
    # plt.tight_layout()
    # sns.despine(fig)
    # plt.savefig('/home/cjensen/Desktop/phospho_proteome_peptide_count_per_batch.png')
    #

    # make one per patient

    # # ----------------------------------------------------------------------------------
    # # violinplot of protein and phospho sites in data
    #
    #
    # df_new = pp_clinically_processed.copy()
    # df_new['Site positions'] = df_new['Site positions'].fillna('nan')
    # psite_list = []
    # for row_index in range(0, df_new.shape[0]):
    #     if psite_list:
    #         psite_list.append(df_new.iloc[row_index,
    #                                       df_new.columns.get_loc(
    #                                           'Site positions')].split(
    #             ';'))
    #     else:
    #         psite_list = df_new.iloc[row_index,
    #                                  df_new.columns.get_loc(
    #                                      'Site positions')].split(';')
    # # psite_list = [item for sublist in psite_list for item in sublist]
    # # print(len(psite_list))
    # # psite_list = set(psite_list)
    # # remove na
    # psite_list = [item for item in psite_list if item != ['nan']]
    # print(len(psite_list))
    #
    #
    # ##  read in Phospho (STY)Sites.txt
    # raw_data_location = '/home/cjensen/kusterlab/internal_projects/active/TOPAS/WP31/Searches/Patients'
    # try:
    #     psites_files = [os.path.join(root, name) for root, dirs, files in
    #                 sorted(os.walk(raw_data_location)) for
    #                 name in files if name.endswith("Phospho (STY)Sites.txt") and
    #                 'PP' in root]
    # except FileNotFoundError:
    #     print("No datafiles found")
    #     sys.exit(2)
    # new_pp_files = psites_files.copy()
    # index = -1
    # for i in range(len(psites_files)):
    #     index += 1
    #     element = psites_files[i]
    #     if 'QC_FAILED' in element:
    #         new_pp_files.pop(index)
    #         index -= 1
    #     if 'outdated' in element:
    #         new_pp_files.pop(index)
    #         index -= 1
    # print(len(new_pp_files))  # 22
    #
    # # sort files after batches
    # sorted_psite_files, batch_numbers = prep.sort_batches(new_pp_files)
    # # print(sorted_psite_files)
    #
    # ## Loop over files
    # psite_list = []
    # for file in sorted_psite_files:
    #     # print(file)
    #     if file == '/home/cjensen/kusterlab/internal_projects/active/TOPAS/WP31' \
    #                '/Searches/Patients/Batch3_PP_MASTER/combined/txt/Phospho (' \
    #                'STY)Sites.txt':
    #         temp = pd.read_csv(file, sep='\t')
    #         # print(temp.columns.tolist())
    #         number_of_phospho = temp.loc[:, 'Number of Phospho (STY)']
    #         print('len with nan', len(number_of_phospho))
    #         # print(number_of_phospho.values.tolist())
    #         # number of lines not being nan
    #         print('len with no nan', len(number_of_phospho.dropna()))
    #         for psite in range(0, number_of_phospho.dropna().shape[0]):
    #             if psite_list:
    #                 psite_list.append(number_of_phospho.dropna().iloc[psite].split(
    #                     ';'))
    #             else:
    #                 # print(psite)
    #                 # print(number_of_phospho.dropna().iloc[psite])
    #                 psite_list = number_of_phospho.dropna().iloc[psite].split(';')
    #         # print(psite_list)
    #         psite_list = [item for sublist in psite_list for item in sublist if item\
    #                                                                               != ['nan']]
    #         psite_int_map = map(int, psite_list)
    #         psite_list = list(psite_int_map)
    #         print('len after split', len(psite_list))
    #         print('sum after split', sum(psite_list))

    # ----------------------------------------------------------------------------------
    # Barplot of protein / phospho peptide numbers in the batches
    #
    # fp_filtered = [135666, 120261, 123858, 94525, 106636, 106592, 80463, 107695, 106607,
    #                100094, 89225, 126064, 124067, 132031, 137316, 115549, 91675, 144221,
    #                147955, 134944, 126995, 132189]
    #
    # pp_filtered = [28235, 26550, 28888, 17720, 24524, 25116, 16800, 22245, 23574, 18826,
    #                15233, 29207, 27607, 26128, 25276, 25759, 34551, 37585, 38285, 29066,
    #                24566, 30598]
    #
    # df = pd.DataFrame({
    #     'Factor': ['Batch 1', 'Batch 2', 'Batch 3', 'Batch 4', 'Batch 5',
    #                'Batch 6', 'Batch 7', 'Batch 8', 'Batch 9', 'Batch 10',
    #                'Batch 11', 'Batch 16', 'Batch 17', 'Batch 18', 'Batch 20',
    #                'Batch 21', 'Batch 22', 'Batch 23', 'Batch 24', 'Batch 25',
    #                'Batch 26', 'Batch 27'],
    #     'Full proteome': fp_filtered,
    #     'Phospho proteome': pp_filtered
    # })
    #
    # fig, ax1 = plt.subplots(figsize=(16, 7))
    # tidy = df.melt(id_vars=['Factor']).rename(columns=str.title)
    # tidy.columns = ['Batch', 'Data', 'Number of peptides']
    # # print(tidy)
    # sns.barplot(x='Batch', y='Number of peptides', hue='Data', data=tidy, ax=ax1)
    # ax1.set_xlabel('Batches', fontsize=14)
    # ax1.set_ylabel('Number of peptides', fontsize=14)
    # ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, fontsize=13)
    # ax1.legend(fontsize=13)
    # # ax1.set_yticklabels(ax1.get_yticklabels(), fontsize=13)
    # plt.tight_layout()
    # sns.despine(fig)
    # plt.savefig('/home/cjensen/Desktop/test.png')

    # Intensity distribution figure showing two patients location

    fgfr1_data = fp_clinically_processed.copy()
    fgfr1_data.drop('gene_names', axis=1, inplace=True)
    # remove reference channels
    if len(fgfr1_data.columns) % 11 == 0:
        n_batches = int(len(fgfr1_data.columns) / 11)
    else:
        print('Error number of columns is {}'.format(len(fgfr1_data.columns)))
    ref_cols = [np.array([8, 9, 10]) + (11 * i) for i in range(n_batches)]
    # flatten list of arrays to list
    ref_cols = np.concatenate(ref_cols).ravel().tolist()
    fgfr1_data.drop(fgfr1_data.columns[ref_cols], axis=1, inplace=True)
    # print(fgfr1_data.index.values)


    # gapdh_data = fgfr1_data.loc['GAPDH', :]
    # # print(gapdh_data)
    abl1_data = fgfr1_data.loc['ABL1', :]
    crkl_data = fgfr1_data.loc['CRKL', :]
    pten_data = fgfr1_data.loc['PTEN', :]
    h12_data = fgfr1_data.loc['H1-2', :]
    fgfr1_data = fgfr1_data.loc['FGFR1', :]

    # # # sort it and plot it
    # df = pd.DataFrame({
    #     'Patients': list(range(0, 176)),
    #     'Normalized intensity': fgfr1_data.sort_values(axis=0, ascending=False)
    # })
    # df2 = pd.DataFrame({
    #     'Patients': list(range(0, 176)),
    #     'Normalized intensity': abl1_data.sort_values(axis=0, ascending=False)
    # })
    # df3 = pd.DataFrame({
    #     'Patients': list(range(0, 176)),
    #     'Normalized intensity': crkl_data.sort_values(axis=0, ascending=False)
    # })
    # df4 = pd.DataFrame({
    #     'Patients': list(range(0, 176)),
    #     'Normalized intensity': pten_data.sort_values(axis=0, ascending=True)
    # })
    # df5 = pd.DataFrame({
    #     'Patients': list(range(0, 176)),
    #     'Normalized intensity': h12_data.sort_values(axis=0, ascending=False)
    # })
    #
    # # print(len(fgfr1_data.values)) # 176
    #
    # with pd.option_context('display.max_rows', None, 'display.max_columns',
    #                        None):  # more options can be specified also
    #     print(df)



    # for case study 1, for fgfr1 use sorted patient 6, for abl1 65, for crkl 1
    # import pylab
    # params = {'legend.fontsize': 'x-large',
    #           'figure.figsize': (5, 5),
    #           'axes.labelsize': 'x-large',
    #           'axes.titlesize': 'x-large',
    #           'xtick.labelsize': 'x-large',
    #           'ytick.labelsize': 'x-large'}
    # pylab.rcParams.update(params)
    # grid = sns.JointGrid(x='Patients', y='Normalized intensity', data=df5)
    # g = grid.plot_joint(sns.scatterplot, data=df5, s=30, marker='o')
    # g.ax_joint.scatter(5, 2.880778, marker='o', s=40, c='r')
    # # # g.ax_joint.scatter(104, 0.002252, marker='o', s=12, c='r')
    # g.ax_joint.set_ylabel('Normalized intensity', fontsize=14)
    # g.ax_joint.set_xlabel('Patients', fontsize=14)
    # g.ax_joint.set_title('H1-2 expression in cohort', fontsize=16)
    # sns.kdeplot(df['Normalized intensity'], ax=g.ax_marg_y, legend=False, vertical=True)
    # g.ax_marg_y.axhline(2.880778, 5, 1, color='r')
    # # # # # median
    # # # # g.ax_marg_y.axhline(0.2241696686575727, 0, 1, color='mediumblue', linewidth=3)
    # # # # g.ax_marg_y.axhline(0.002252, 0, 1, color='r')
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/rank_plot_h12.png')

    # sns.scatterplot(x='Patients', y='Normalized intensity', data=df4, s=20, marker='o')
    # plt.scatter(x=36, y=0.030025, marker='o', s=40, c='r')
    # plt.title('PTEN expression in cohort', fontsize=16)
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/rank_plot_pten2.png')
    #



    # grid = sns.JointGrid(x='Patients', y='Normalized intensity', data=df)
    # g = grid.plot_joint(sns.scatterplot, data=df, s=12, marker='o')
    # g.ax_joint.scatter(0, 1.1771449280794766, marker='o', s=12, c='r')
    # # g.ax_joint.scatter(104, 0.002252, marker='o', s=12, c='r')
    # g.ax_joint.set_ylabel('Normalized intensity', fontsize=14)
    # g.ax_joint.set_xlabel('Patients', fontsize=14)
    # g.ax_joint.set_title('FGFR1 expression in cohort', fontsize=16)
    # sns.kdeplot(df['Normalized intensity'], ax=g.ax_marg_y, legend=False, vertical=True)
    # g.ax_marg_y.axhline(1.1771449280794766, 0, 1, color='r')
    # # median
    # g.ax_marg_y.axhline(0.2241696686575727, 0, 1, color='mediumblue', linewidth=3)
    # # g.ax_marg_y.axhline(0.002252, 0, 1, color='r')
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/rank_plot_fgfr1_1.png')
    # plt.closefig()
    # example patients
    # print(max(df['Normalized intensity']))  # 1.1771449280794766
    # print(df['Normalized intensity'].median())  # 0.2241696686575727
    # print(df['Normalized intensity'].mean())  # 0.2241696686575727
    #
    # print(fgfr1_data[fgfr1_data == 1.1771449280794766])  # Reporter intensity
    # # corrected 5 Batch21_FP_MASTER
    # # print(fgfr1_data[fgfr1_data == 0.2241696686575727]) #
    # # print(fgfr1_data[fgfr1_data.round(1) == 0.002252])
    # print(df[df.index == 'Reporter intensity corrected 6 Batch27_FP_MASTER'])
    # # Reporter intensity corrected 6 Batch27_FP_MASTER

    # sort it and plot it
    # df = pd.DataFrame({
    #     'Patients': list(range(0, 176)),
    #     'Normalized intensity': gapdh_data.sort_values(axis=0, ascending=False)
    # })
    # with pd.option_context('display.max_rows', None, 'display.max_columns',
    #                        None):  # more options can be specified also
    #     print(df)
    # print(len(fgfr1_data.values)) # 176

    # grid = sns.JointGrid(x='Patients', y='Normalized intensity', data=df)
    # g = grid.plot_joint(sns.scatterplot, data=df, s=12, marker='o')
    # g.ax_joint.scatter(0, 0.8093648286100272, marker='o', s=12, c='r')
    # g.ax_joint.scatter(73, -0.015825, marker='o', s=12, c='r')
    # g.ax_joint.set_ylabel('Normalized intensity', fontsize=14)
    # g.ax_joint.set_ylim((-2, 1.5))
    # g.ax_joint.set_xlabel('Patients', fontsize=14)
    # g.ax_joint.set_title('GAPDH expression in cohort', fontsize=16)
    # sns.kdeplot(df['Normalized intensity'], ax=g.ax_marg_y, legend=False,
    #             vertical=True)
    # g.ax_marg_y.axhline(0.8093648286100272, 0, 1, color='r')
    # # median
    # g.ax_marg_y.axhline(-0.049652874261551155, 0, 1, color='mediumblue',
    #                     linewidth=3)
    # # g.ax_marg_y.axhline(0.002252, 0, 1, color='r')
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/rank_plot_gapdh.png')
    # example patients
    # print(max(df['Normalized intensity']))  # 0.8093648286100272
    # print(df['Normalized intensity'].median())  # -0.049652874261551155
    # print(df['Normalized intensity'].mean())  # -0.057781270275141995

    # Reporter intensity corrected 8 Batch10_FP_MASTER
    # Reporter intensity corrected 4 Batch27_FP_MASTER
    # cols = ['darkblue', '#ADD8E6', ]
    # labels = ['identified', 'non-identified']
    # patches, texts = plt.pie([337, 540 - 337], explode=(0.05, 0.05),
    #         startangle=90, radius=0.8,
    #         colors=cols)
    # plt.legend(patches, labels, loc='lower center', ncol=2, frameon=False, fontsize=16)
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/pie_chart_oncogenes_ident.png')
    #
    # patches, texts = plt.pie([297, 540 - 297], explode=(0.05, 0.05),
    #         startangle=90, radius=0.8,
    #         colors=cols)
    # plt.legend(patches, labels, loc='lower center', ncol=2, frameon=False, fontsize=16)
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/pie_chart_oncogenes_ident_psite.png')



    # # Read in and plot the p-peptides vs protein numbers per sample

    import pylab
    params = {'legend.fontsize': 'x-large',
              'figure.figsize': (7, 5),
              'axes.labelsize': 'x-large',
              'axes.titlesize': 'x-large',
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'}
    pylab.rcParams.update(params)
    df = pd.read_csv('/home/cjensen/kusterlab_home/kusterlab/file_exchange/forCecilia/fromJulia'
                '/protein_peptide_numbers.csv')
    print(df)
    df.rename({'sites': 'phospho-peptides'}, axis=1, inplace=True)
    print(df)
    test = df[['proteins', 'phospho-peptides']].stack().reset_index()

    df2 = pd.DataFrame({
        'Category': test.loc[:, 'level_1'].values,
        'Numbers': test[0].values
    })

    g = sns.swarmplot(x='Category', y='Numbers', data=df2)
    plt.ylim(0, 27000)
    ylabels = [str(int(y)) + 'K'  for y in g.get_yticks() / 1000
               ]
    g.set_yticklabels(ylabels)
    plt.tight_layout()
    # plt.xlabel(fontsize=14)
    plt.savefig('/home/cjensen/Desktop/protein_peptide_swarmplot3.png')


    # plot for patient 1 ranks  - makes no sense
    #
    # rank_case1 = pd.read_csv('/home/cjensen/kusterlab_home/kusterlab/internal_projects/active/TOPAS'
    #                          '/WP31/Analyses/MTB_reports/2021.10.26/full_proteome_measures.csv')
    #
    # # print(rank_case1.columns.values[1:-1])
    # # print(rank_case1[rank_case1['Protein IDs'] == 'FGFR1'].values[0][1:-1])
    # names = rank_case1.columns.values[1:-1]
    # values = rank_case1[rank_case1['Protein IDs'] == 'FGFR1'].values[0][1:-1]
    #
    # # sort values and sort names with sorting index
    #
    # df = pd.DataFrame({
    #     'Patients labels': names,
    #     'Normalized intensity': values
    # })
    # print(df)
    # df = df.sort_values(by='Normalized intensity', axis=0, ascending=False)
    # print(df)
    #
    # with pd.option_context('display.max_rows', None, 'display.max_columns',
    #                        None):  # more options can be specified also
    #     print(df)

    # grid = sns.JointGrid(x='Patients', y='Normalized intensity', data=df)
    # g = grid.plot_joint(sns.scatterplot, data=df, s=12, marker='o')
    # g.ax_joint.scatter(0, 1.1771449280794766, marker='o', s=12, c='r')
    # g.ax_joint.scatter(104, 0.002252, marker='o', s=12, c='r')
    # g.ax_joint.set_ylabel('Normalized intensity', fontsize=14)
    # g.ax_joint.set_xlabel('Patients', fontsize=14)
    # g.ax_joint.set_title('FGFR1 expression in cohort', fontsize=16)
    # sns.kdeplot(df['Normalized intensity'], ax=g.ax_marg_y, legend=False, vertical=True)
    # g.ax_marg_y.axhline(1.1771449280794766, 0, 1, color='r')
    # # median
    # g.ax_marg_y.axhline(0.2241696686575727, 0, 1, color='mediumblue', linewidth=3)
    # g.ax_marg_y.axhline(0.002252, 0, 1, color='r')
    # plt.tight_layout()
    # plt.savefig('/home/cjensen/Desktop/rank_plot.png')