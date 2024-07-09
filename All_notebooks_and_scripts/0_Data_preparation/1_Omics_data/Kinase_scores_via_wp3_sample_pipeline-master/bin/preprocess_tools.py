import re
import os
from typing import List, Dict, Tuple, Union, Callable
from pathlib import Path
import logging

import pandas as pd
import numpy as np

from . import utils
from . import sample_metadata
from . import sample_annotation
from . import identification_metadata as id_meta

logger = logging.getLogger(__name__)
MQ_EVIDENCE_COLUMNS = ['Modifications', 'Modified sequence',
                       'Proteins', 'Leading proteins', 'Gene names',
                       'Type', 'Raw file', 'Fraction', 'Experiment',
                       'Charge', 'PEP', 'MS/MS scan number', 'Score', 'Intensity',
                       'Reverse', 'Potential contaminant', 'id'] + \
                      [f'Reporter intensity corrected {i}' for i in range(1, 12)]

MQ_EVIDENCE_COLUMNS_TYPES = {'Modifications': 'category', 'Modified sequence': '', 'Proteins': '', 'Leading proteins': '',
                             'Gene names': '', 'Type': '', 'Raw file': '', 'Fraction': '', 'Experiment': 'category', 'Charge': '',
                             'PEP': '', 'MS/MS scan number': '', 'Score': '', 'Intensity': '',
                             'Reverse': '', 'Potential contaminant': '', 'id': ''}


def check_annot(sample_annotation_file: str, metadata_annotation_file: str, in_metadata: Callable, patient_regex: str):
    """Performs consistency checks of the sample annotation and metadata files.

    Args:
        sample_annotation_file: _description_
        metadata_annotation_file: _description_
        in_metadata: _description_
        patient_regex: _description_

    Raises:
        ValueError: _description_
        ValueError: _description_
        ValueError: _description_
        ValueError: _description_

    Returns:
        _description_
    """    
    sample_annot_df = sample_annotation.load_sample_annotation(sample_annotation_file)
    sample_annot_df_filtered = sample_annotation.filter_sample_annotation(sample_annot_df, remove_qc_failed=True, remove_replicates=True)
    sample_annot_df_filtered = sample_annot_df_filtered.reset_index(drop=True)
    metadata_df = sample_metadata.load(metadata_annotation_file)
    # TODO: get patients not matching regex and print in error

    # Check if all patient annot fit with patient regex
    if sample_annot_df_filtered.set_index(['Sample name']).filter(regex=patient_regex, axis=0).shape[0] != sample_annot_df_filtered['Sample name'].shape[0]:
        raise ValueError('Not all samples match given patient regex')

    # Check that there are no duplicates in neither patient annot and metadata
    if sample_annot_df_filtered['Sample name'].duplicated().any():
        duplicated = sample_annot_df_filtered[sample_annot_df_filtered['Sample name'].duplicated()]
        raise ValueError(f'Duplicated sample(s) in sample annotation: {duplicated}')

    elif metadata_df['Sample name'].duplicated().any():
        duplicates_samples =  metadata_df['Sample name'][metadata_df['Sample name'].duplicated()]
        duplicates_samples = ('_').join(list(duplicates_samples))
        raise ValueError(f'The Duplicated sample(s): {duplicates_samples} in metadata annotation should be removed')

    # TODO : needs to be refactored

    # Check if all given samples in patient annot is also in metadata
    if sample_annot_df_filtered['Sample name'][~sample_annot_df_filtered['Sample name'].apply(lambda x: in_metadata(x)).isin(metadata_df['Sample name'])].any():
        missing_samples = sample_annot_df_filtered['Sample name'][~sample_annot_df_filtered['Sample name'].apply(lambda x: in_metadata(x)).isin(metadata_df['Sample name'])].values
        #raise ValueError(f'Sample(s) not found in metadata: {missing_samples}')
        print(f'Sample(s) not found in metadata: {missing_samples}')

    return sample_annot_df
    

def in_metadata(sample: str = 'CHD'):
    if sample.split('-')[-1] == 'R2':
        return '-'.join(sample.split('-')[:-1])
    return sample


def remove_ref_empty_batch(df_with_ref: pd.DataFrame, sample_annotation_df: pd.DataFrame):
    # remove ref for empty batches
    sample_names_drop = []
    for col in df_with_ref.columns:
        if 'Reporter intensity corrected' in col:
            batch = re.search('Batch\d{1,2}', col).group()
            batch = re.findall(r'\d+', batch)[0]
            qc_passed = sample_annotation_df.loc[sample_annotation_df['Batch Name'] == int(batch), 'QC']
            if 'passed' not in qc_passed.values and 'shaky' not in qc_passed.values:
                sample_names_drop.append(col)
    # samples have already been dropped but we want to drop the reporter intensity ones for this batch
    df_with_ref = df_with_ref.drop(columns=sample_names_drop)
    return df_with_ref


def get_evidence_files(sample_annotation_df: pd.DataFrame,
                       raw_data_location: Union[str, Path],
                       data_type: str) -> List:
    batches = sample_annotation.get_unique_batches(sample_annotation_df)

    evidence_files = get_data_location(raw_data_location, data_type)

    # remove path if batch not in annotation file
    evidence_files = filter_evidence_files(evidence_files, data_type.upper(), batches)
    return evidence_files


def log_transform_intensities(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Log10 transforming intensities")
    tmt_channels = utils.get_tmt_channels(df)
    df.loc[:, tmt_channels.columns] = np.log10(tmt_channels.replace(0, np.nan))
    return df


def sum_peptide_intensities(df: pd.DataFrame,
                            debug: bool = False, run_lfq: bool = False) -> pd.DataFrame:
    """sum up intensities per p-peptide across fractions and charge states"""
    logger.info("Summing up intensities per p-peptide across fractions and charge states")
    # TODO: split strings on semicolon before making unique
    def csv_list_unique(x):
        return ";".join(map(str, list(dict.fromkeys(x))))

    def aggregate_imputations(x):
        annotations = dict.fromkeys(x)
        if not set(annotations).issubset({'imputed;', ''}):
            raise ValueError(
                f"Found other annotations ({annotations}) besides imputations, need new solution to detect partially imputed peptides")

        if set(annotations) == {'imputed;', ''}:
            return 'partially imputed;'

        return ";".join(map(str, list(annotations)))

    # Phospho preprocessing
    # if 'Modified sequence' in df.columns and debug and run_lfq is False:
    #     print('hello')
    #     print(df.columns)
    #     df = df.groupby(['Batch', 'Modified sequence', 'Gene names'])
    #     print('step 1')
    #
    #     df = df.agg(
    #         # **{
    #         #     'new_genes': pd.NamedAgg(column='new_genes', aggfunc=csv_list_unique),
    #         # },
    #         **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc=sum) for i in
    #            range(1, 12)}
    #     ).reset_index()
    #     print('step 2')
    # Phospho raw for last module? TODO check this
    if 'Modified sequence' in df.columns and run_lfq is False:
        df = df.groupby(['Batch', 'Modified sequence'])

        # if 'Transferred spectra count' in df.columns:
        #     df = df.agg(
        #         **{
        #             'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique),
        #             'Gene names': pd.NamedAgg(column='Gene names', aggfunc=csv_list_unique),
        #             'Transferred spectra count': pd.NamedAgg(column='Transferred spectra count', aggfunc=csv_list_unique),
        #         },
        #         **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc=sum) for i
        #            in
        #            range(1, 12)},
        #         **{f'Identification metadata {i}': pd.NamedAgg(column=f'Identification metadata {i}',
        #                                                        aggfunc=aggregate_imputations)
        #            for i in
        #            range(1, 12)}
        #     ).reset_index()
        # else:
        df = df.agg(
            **{
                'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique),
                'Gene names': pd.NamedAgg(column='Gene names', aggfunc=csv_list_unique),
            },
            **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc=sum) for i in
               range(1, 12)},
            **{f'Identification metadata {i}': pd.NamedAgg(column=f'Identification metadata {i}', aggfunc=aggregate_imputations)
               for i in
               range(1, 12)}
        ).reset_index()
    elif 'Modified sequence' in df.columns and run_lfq:
        df = df.groupby(['Batch', 'Modified sequence'])
        df = df.agg(
            **{
                'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique),
                'Gene names': pd.NamedAgg(column='Gene names', aggfunc=csv_list_unique),
            },
            **{'Reporter intensity corrected 1': pd.NamedAgg(column='Reporter intensity corrected 1', aggfunc=sum)},
            **{'Identification metadata 1': pd.NamedAgg(column='Identification metadata 1', aggfunc=aggregate_imputations)}
        ).reset_index()
    else:
        df = df.groupby(['Batch', 'Gene names'])
        df = df.agg(
            **{
                'Proteins': pd.NamedAgg(column='Proteins', aggfunc=csv_list_unique),
                'new_genes': pd.NamedAgg(column='new_genes', aggfunc=csv_list_unique),
            },
            **{f'Reporter intensity corrected {i}': pd.NamedAgg(column=f'Reporter intensity corrected {i}', aggfunc=sum) for i in
               range(1, 12)}
        ).reset_index()
    return df


def load_and_normalize(data_loader: 'DataLoader',
                       results_folder: Union[str, Path],
                       data_type: str,
                       debug: bool = True) -> pd.DataFrame:
    """
    Function for preprocessing MQ evidence files
    :param data_loader:
    :param results_folder: location of log and intermediately processed files
    :param data_type:
    :param debug: save extra output useful for debugging
    :return: Dataframe with combined preprocessed data from all batches
    """
    # Parse evidence files
    all_batches = data_loader.load_data(MQ_EVIDENCE_COLUMNS)

    if debug:
        df_raw = pd.concat(all_batches, ignore_index=True)
        df_raw.to_csv(os.path.join(results_folder, f'debug_preprocessed_{data_type}_complete_raw.csv'), index=False)

    # Inside batch median centering
    all_batches = data_loader.median_centering_within_batch(all_batches)
    if debug:
        df_after_first_median = pd.concat(all_batches, ignore_index=True)
        df_after_first_median.to_csv(os.path.join(results_folder, f'debug_preprocessed_{data_type}_after_1st_median.csv'),
                                     index=False)
    # MS1 median centering
    all_batches = data_loader.median_centering_ms1(all_batches)

    # Combine all batches into a single dataframe
    df = pd.concat(all_batches, ignore_index=True)
    if debug:
        df.to_csv(os.path.join(results_folder, f'debug_preprocessed_{data_type}_after_ms1_centering.csv'), index=False)

    # Impute an MS1 value for scans without MS1 based on mean MS1 share of reference channels in other batches
    df = data_loader.impute_ms1_intensity(df)
    if debug:
        df.to_csv(os.path.join(results_folder, f'debug_preprocessed_{data_type}_after_ms1_imputation.csv'), index=False)

    # Distribute MS1 intensity over the TMT channels
    df = data_loader.redistribute_ms1_intensity(df)
    if debug:
        df.to_csv(os.path.join(results_folder, f'debug_preprocessed_{data_type}_after_ms1_correction.csv'), index=False)

    return df


def convert_long_to_wide_format(df: pd.DataFrame,
                                has_metadata_cols: bool = False) -> pd.DataFrame:
    """
    Converts dataframe from long format (11 reporter channels, batches concatenated row-wise)
    to wide format (samples as columns)
    """
    logger.info("Converting from long to wide format")

    tmt_channels = utils.get_tmt_channels(df).columns
    keep_columns = tmt_channels.tolist()

    if has_metadata_cols:
        metadata_cols = id_meta.as_metadata_columns(tmt_channels.str)
        keep_columns += metadata_cols.tolist()

    if 'Modified sequence' in df.columns:
        keep_columns += ['Batch', 'Gene names', 'Proteins']
        if 'Transferred spectra count' in df.columns:
            keep_columns += ['Transferred spectra count']

        index_col = 'Modified sequence'
    else:
        keep_columns += ['Batch']
        index_col = ['Gene names', 'Proteins']
    indexed_df = df.set_index(index_col)[keep_columns]

    all_batches = list()
    for batch_name, df in indexed_df.groupby('Batch'):
        df = df.drop(columns='Batch')
        df = df.rename(columns=lambda x: f'{x} {batch_name}')
        all_batches.append(df)
    wide_df = pd.DataFrame().join(all_batches, how="outer")
    wide_df = wide_df.reset_index()

    if has_metadata_cols:
        wide_df = aggregate_csv_columns(wide_df, 'Gene names')
        wide_df = aggregate_csv_columns(wide_df, 'Proteins')
    return wide_df


def aggregate_csv_columns(df: pd.DataFrame,
                          column_name: str) -> pd.DataFrame:
    all_cols = df.columns[df.columns.str.startswith(column_name)]

    def csv_list_unique(csv_string):
        l = csv_string.split(';')
        s = sorted({x for x in l if x != 'nan'})
        return ";".join(s)

    df[column_name] = df[all_cols].astype(str).agg(';'.join, axis=1).apply(csv_list_unique)

    df = df.drop(columns=all_cols)
    return df


def drop_duplicate_indices(df: pd.DataFrame) -> pd.DataFrame:
    """Keeps only first entry for entries with duplicated indices"""
    return df[~df.index.duplicated(keep='first')]


def merge_ms1_columns(all_batches: List[pd.DataFrame]) -> pd.DataFrame:
    ms1s = list()
    for df in all_batches:
        batch_name = df['Batch'].iloc[0]
        ms1_intensity_df = drop_duplicate_indices(
            df.set_index('Modified sequence')[['Intensity']].rename(columns={'Intensity': batch_name}))
        ms1s.append(ms1_intensity_df)

    # join MS1 intensities by modified sequence
    return pd.DataFrame().join(ms1s, how="outer")


def median_centering(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes samples by multiplying each sample with its own correction factor
    calculated by division of average median of reference channels with the median of
    the sample. That way sample distributions gets centered around the same value
    :param df: dataframe of batch intensities (11 reporter channels) or dict of MS1 values per batch
    """
    # Multiplying with correction factor to center samples
    df = df.replace(0, np.nan)
    medians = df[df.isna().sum(axis=1) <= 3].median(axis=0)
    df = df * (medians.mean() / medians)
    return df


def impute_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    If a peptide is not detected for a channel but is detected in the batch,
    impute with 1/100th of the maximum intensity or the minimum intensity for 
    a peptide within its batch, whichever is lower.
    
    :param df: dataframe of batch intensities (11 reporter channels)
    """
    logger.info('Imputing data')

    tmt_cols_df = utils.get_tmt_channels(df)
    tmt_cols_df = tmt_cols_df.replace(0, np.nan)
    patient_channels_df = tmt_cols_df.filter(regex=r'^Reporter intensity corrected [1-8]$')

    x = tmt_cols_df.max(axis=1) / 100
    y = tmt_cols_df.min(axis=1)
    imputation_values = pd.concat([x, y], axis='columns').min(axis='columns')

    # apply per column
    patient_channels = patient_channels_df.columns
    metadata_cols = id_meta.as_metadata_columns(patient_channels.str)

    def add_imputation_status(z):
        return np.where(z.isna() & imputation_values.notna(), 'imputed;', '')

    df[metadata_cols] += tmt_cols_df.apply(add_imputation_status).rename(columns=id_meta.as_metadata_columns)
    df.loc[:, patient_channels] = patient_channels_df.apply(lambda z: z.fillna(imputation_values))
    return df


def filter_data(df: pd.DataFrame,
                data_type: str) -> pd.DataFrame:
    """
    filters out contaminants, reverse sequences from decoy database,
    For phospho, non-phospho modified peptides are removed as well.
    :param df:
    :param results_folder:
    :param data_type:
    """
    logger.info(f'Filtering data {data_type}')

    logger.info(f'- before filtering: {df.shape[0]}')

    df = df[df['Potential contaminant'] != '+']
    logger.info(f'- after contaminant removal: {df.shape[0]}')

    df = df[df['Reverse'] != '+']
    logger.info(f'- after reverse removal: {df.shape[0]}')

    df = df.drop(['Reverse', 'Potential contaminant'], axis=1)

    if data_type == 'pp':
        df = df[df['Modifications'].str.contains('Phospho (STY)', regex=False)]
        logger.info(f'- after unmodified peptides removal: {df.shape[0]}')

        df = df.drop('Modifications', axis=1)

    return df


def rename_columns_with_sample_ids(df: pd.DataFrame,
                                   channel_to_sample_id_dict: Dict[str, str],
                                   index_cols: List[str],
                                   remove_ref=True) -> pd.DataFrame:
    """
    Transform column names of the format 'Reporter intensity corrected <TMT_channel> <batch_name>' to the sample names
    """
    df = df.set_index(index_cols)

    tmt_channels = list(channel_to_sample_id_dict.keys())
    metadata_cols = list(map(id_meta.as_metadata_columns, tmt_channels))

    keep_cols = tmt_channels + metadata_cols

    if not remove_ref:
        ref_channels_cols = df.filter(regex='^Reporter intensity corrected (9|10|11)').columns.tolist()
        keep_cols += ref_channels_cols

    df = df.filter(items=keep_cols, axis=1)

    # build dictionary to also rename the metadata columns with the sample ids
    metadata_cols_with_sample_id = map(lambda x: f'Identification metadata {x}', channel_to_sample_id_dict.values())
    metadata_to_sample_id_dict = dict(zip(metadata_cols, metadata_cols_with_sample_id))

    rename_dict = {**channel_to_sample_id_dict, **metadata_to_sample_id_dict}
    df = df.rename(columns=rename_dict)
    df = df.reset_index()
    return df


def get_data_location(maxquant_super_folder: Union[str, Path],
                      data_type: str,
                      file_type="evidence.txt") -> Tuple[List[str], List[str]]:
    """
    Return file paths of all relevant files from maxquant_super_folder folder
    :param maxquant_super_folder: Folder containing results from MQ searches
    :param data_type:
    :param file_type:
    :return: List of file paths for full and phospho proteome data
    """
    if not os.path.isdir(maxquant_super_folder):
        raise IOError("No datafiles found")

    validity_check = is_valid_fp_file
    if data_type == 'pp':
        validity_check = is_valid_pp_file

    evidence_files = []
    for directory, _, files in sorted(os.walk(maxquant_super_folder)):
        for filename in files:
            if not filename.endswith(file_type):
                continue

            path = os.path.join(directory, filename)
            if validity_check(directory):
                evidence_files.append(path)
    return evidence_files


def filter_evidence_files(files: List[Union[str, Path]],
                          data_type: str,
                          batches: List[str]):

    # test that we can crash it if multiple batch 40

    return [file for file in files for batch_number, cohort in batches if
            os.path.join(cohort, f'Batch{batch_number}_{data_type}') in file]


def is_valid_fp_file(directory: Union[str, Path]) -> bool:
    return 'FP' in directory \
        and 'Merged_FP' not in directory \
        and 'QC_FAILED' not in directory \
        and 'outdated' not in directory \
        and 'Technician' not in directory \
        and '_var_' not in directory


def is_valid_pp_file(directory: Union[str, Path]) -> bool:
    return 'PP' in directory \
        and 'QC_FAILED' not in directory \
        and 'outdated' not in directory \
        and 'Technician' not in directory \
        and '_var_' not in directory