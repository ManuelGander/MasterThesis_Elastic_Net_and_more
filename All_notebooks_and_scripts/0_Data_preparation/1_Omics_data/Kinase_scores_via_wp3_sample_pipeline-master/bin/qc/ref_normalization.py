import os
import re
import sys
import argparse
import logging
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import bin.config as config

logger = logging.getLogger(__name__)


def perform_normalization(df: pd.DataFrame, batches: List, ref: List):
    """ Code translated from R and written by Chen"""
    # ref = None

    # TODO: test with ref channels
    # TODO: check if missing values introduced
    # TODO: if missing values introduced, impute missing ref channel values


    if ref is None:
        ref = np.arange(len(df.columns))
    # batches_ref = [batches[i] for i in ref]
    batches_ref = batches[ref]
    # print(batches_ref)
    df_ref = df.loc[:, np.array(ref)]
    # print(df_ref)

    overall_means = df.median(axis=1)
    for batch in batches.unique():
        # print(batch)  # 1
        batch_ref_bool = batches_ref == batch
        batch_bool = batches == batch
        # print(batch_ref_bool)  # 177
        # batch_samples = df.loc[:, np.array(batch_ref_bool)].columns.tolist()
        batch_samples = df_ref.loc[:, np.array(batch_ref_bool)].columns.tolist()
        # print(batch_samples)
        off_set = df_ref.loc[:, np.array(batch_ref_bool)].median(axis=1) - overall_means
        # print(off_set)
        df.iloc[:, np.array(batch_bool)] = df.iloc[:, np.array(batch_bool)].sub(off_set, axis=0, level='Modified sequence')
    # print(df.shape)
    return df, batches_ref


def read_data(results_folder: str, sample_annotation: str):
    df = pd.read_csv(os.path.join(results_folder, 'preprocessed_pp_with_ref.csv'))
    sample_annot = pd.read_csv(sample_annotation)
    df = df.set_index(['Modified sequence', 'Gene names', 'Proteins'])
    # filter out ident metadata
    # df = df.loc[:, ~df.columns.str.startswith('Identification')]
    return df, sample_annot


def normalize_df(df: pd.DataFrame, sample_annot: pd.DataFrame):
    # create batch list
    batches, refs = create_batch_number_list(df, sample_annot)
    # new_df = pd.melt(df, ignore_index=False)
    # zipped_pairs = dict(zip(df.columns.tolist(), batches))
    # new_df['Batch'] = new_df['variable'].apply(lambda x: zipped_pairs[x])

    # sort list and df accordingly - channel order not 9,10,11 but should not matter - order is 10, 11, 9
    sorted_df, batches, refs = sort_after_batches(df, batches, refs)
    # print(sorted_df.columns)
    # print(batches)
    # print(refs)
    # print(df.filter(df[df.columns.str.startswith('Reporter intensity')].columns, axis=1).describe())

    # sum per batch of missing ref channels
    missingness_df = pd.DataFrame()
    for i in range(0, len(sorted_df.loc[:, np.array(refs)].columns), 3):
        bool_df = sorted_df.loc[:, np.array(refs)].isna()
        # if 2 in bool_df.iloc[:, i:i+3].sum(axis=1):
            # print(bool_df.iloc[:, i:i+3].sum(axis=1))
        if i == 0:
            missingness_df = bool_df.iloc[:, i:i+3].sum(axis=1)
        else:
            missingness_df = pd.concat([missingness_df, bool_df.iloc[:, i:i+3].sum(axis=1)], axis=1)
    # print(missingness_df)
    missingness_df.to_csv('/home/cjensen/Desktop/missingness_df.csv')
    x = missingness_df.sum(axis=1)/3
    fig = plt.figure()
    ax = x.hist(bins=59)
    ax.set_xlabel("number of batches")
    ax.set_ylabel("missingness of p-peptides in ref channels")
    fig.savefig('/home/cjensen/Desktop/histogram_missingness_ref_channels_pp.png', dpi=600)
    #
    # print(sorted_df.shape)

    # sorted_filtered_df = sorted_df.drop(sorted_df[missingness_df.sum(axis=1)/3 == 59].index)
    # df = df.drop(sorted_df[missingness_df.sum(axis=1)/3 == 59].index)

    # check for need of performing normalization
    # count_all, count_three_ref_missing = 0, 0
    # for batch_name in batches.unique():
    #     if batch_name == 1:
    #         print('Processing Batch #', batch_name)
    #
    #         # get batch
    #         temp_bool = np.array(batches == batch_name)
    #         temp_batch = sorted_filtered_df.loc[:, temp_bool]
    #         print(temp_batch.iloc[:, -3:].isna().sum(axis=1))
    #         print(temp_batch.iloc[:, -3:].isna().sum(axis=1) == 3)
    #         print(pd.Series(temp_batch.iloc[:, -3:].isna().sum(axis=1) == 3).sum())

            # temp_batch.apply(lambda x: print(x[-3:]), axis=1)
    # print(sorted_filtered_df.shape)
    # print(sorted_filtered_df.isna().sum(axis=1).sum())
    sorted_filtered_normalized_df, batches = perform_normalization(sorted_df, batches, refs)
    # print(df.head())

    # missingness introduced
    # print(df.isna().sum(axis=1).sum())

    # missingness_df = pd.DataFrame()
    # for i in range(0, len(df.loc[:, np.array(refs)].columns), 3):
    #     bool_df = df.loc[:, np.array(refs)].isna()
    #     if i == 0:
    #         missingness_df = bool_df.iloc[:, i:i + 3].sum(axis=1)
    #     else:
    #         missingness_df = pd.concat([missingness_df, bool_df.iloc[:, i:i + 3].sum(axis=1)], axis=1)
    #
    # y = missingness_df.sum(axis=1)/3
    # fig = plt.figure()
    # y.hist(bins=59)
    # fig.savefig('/home/cjensen/Desktop/histogram_missingness_ref_channels_pp_after.png', dpi=600)
    # # simple boxplot
    # batches = batches[:55]
    # df = df.iloc[:, :55]
    # df['Batch'] = pd.Series([np.repeat(batch, df.shape[0]) for batch in batches])
    # fig = plt.figure()
    # boxplot = df.boxplot()
    # fig.savefig('/home/cjensen/Desktop/boxplot_row_wise_normalization.png', dpi=600)

    # sort columns according to order before
    sorted_filtered_normalized_df = sorted_filtered_normalized_df.reindex(columns=df.loc[:, ~df.columns.str.startswith('Identification')].columns)
    # add metadata information again
    df = pd.concat([sorted_filtered_normalized_df, df.loc[:, df.columns.str.startswith('Identification')]], axis=1)
    return df


def create_batch_number_list(df: pd.DataFrame, sample_annot: pd.DataFrame):
    
    # filter out ident metadata
    temp_df = df.loc[:, ~df.columns.str.startswith('Identification')]
    
    batches = pd.Series(temp_df.columns.tolist())
    refs = pd.Series(temp_df.columns.tolist())

    for i, sample in enumerate(batches):
        if sample in sample_annot['Sample name'].values:
            batch_index = sample_annot[sample_annot['Sample name'] == sample].index.values
            batches[i] = sample_annot.loc[batch_index[0], 'Batch Name']
            refs[i] = False
        else:
            if 'Reporter intensity corrected' in sample:
                batch = re.search('Batch\d{1,2}', sample).group()
                batch = re.findall(r'\d+', batch)[0]
                batches[i] = int(batch)
                refs[i] = True
            else:
                print('error')
    return batches, refs


def sort_after_batches(df: pd.DataFrame, batches, refs):
    
    # filter out ident metadata
    temp_df = df.loc[:, ~df.columns.str.startswith('Identification')]
    
    zipped_pairs = zip(batches, refs, temp_df.columns.tolist())
    batches, refs, temp_df_columns = list(zip(*[(batch, ref, column) for batch, ref, column in sorted(zipped_pairs)]))
    temp_df = temp_df.reindex(columns=temp_df_columns)
    return temp_df, pd.Series(batches), pd.Series(refs)


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", dest="results", default=None, metavar="DIR", required=True, help="Path to results to check")

    argv = sys.argv[1:]
    args = parser.parse_args(argv)
    results_folder = Path(args.results)

    config_file = os.path.join(results_folder, 'configs.json')
    configs = config.load(config_file)

    df, sample_annot = read_data(results_folder, configs['sample_annotation'])
    # print(df.columns)

    df = normalize_df(df, sample_annot)


    # second type of normalization?


    # add metadata information again
    # add metadata information again

    # save normalized data with ref channels


    # # remove ref channels
    df = df.drop(df.loc[:, df.columns.str.contains('Reporter intensity')].columns, axis=1)

    # save normalized data
    df.to_csv(os.path.join(results_folder, 'normalized_preprocessed_pp.csv'))

