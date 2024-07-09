import os
import sys
import logging

import pandas as pd
from typing import Union, List
from pathlib import Path

import bin.preprocess_tools as prep
import bin.picked_group as picked
from . import sample_annotation
from . import utils
from . import identification_metadata as id_meta
from bin.data_loaders.tmt_loader import TMTLoader
from bin.data_loaders.simsi_tmt_loader import SimsiTMTLoader
from bin.data_loaders.lfq_loader import LFQLoader

logger = logging.getLogger(__name__)


def preprocess_raw(*args, **kwargs) -> None:
    data_types = kwargs.pop('data_types')
    for data_type in data_types:
        preprocess_raw_data_type(*args, **kwargs, data_type=data_type)


def preprocess_raw_data_type(results_folder: Union[str, Path],
                             sample_annotation_file: Union[str, Path],
                             metadata_annotation: Union[str, Path],
                             patient_regex: str,
                             simsi_folder: Union[str, Path],
                             raw_data_location: Union[str, Path],
                             picked_fdr: float,
                             fasta_file: str,
                             fdr_num_threads: int,
                             program: str,
                             entity: str,
                             histologic_subtype: str,
                             imputation: bool,
                             run_simsi: bool,
                             run_lfq: bool,
                             debug: bool,
                             data_type: str) -> None:
    """
    Calling function for preprocessing of both phospho and full proteome data
    :param results_folder: location in which to save preprocessed files
    :param sample_annotation_file: file containing possible subfolder (of raw folder), metadata and sample information (tmt: batch and channel)
    :param simsi_folder: location in which simsi saves its processed files
    :param raw_data_location: location of data to preprocess (MaxQuant search folders)
    :param picked_fdr: FDR value for picked group filtering
    :param program: subsetting of data for <program> during rest of pipeline when multiple is used for simsi and normalization
    :param entity: subsetting of data for <entity> during rest of pipeline when multiple is used for simsi and normalization
    :param histologic_subtype: subsetting of data for <histologic subtype> during rest of pipeline when multiple is used for simsi and normalization
    :param imputation: boolean for whether to impute or not (only for TMT phospho)
    :param run_simsi: boolean for whether to use simsi for identification transfer between batches
    :param run_lfq: boolean for when data is lfq else tmt is expected
    :param debug: boolean for saving more intermediate results files for debugging purpose
    :param data_type: full protome (fp) or phospho proteome (pp) for which data is analyzed (for lfq only fp)
    """

    # check if file is there - if so skip this
    if os.path.exists(os.path.join(results_folder, f'preprocessed_{data_type}.csv')):
        logger.info(f'Preprocessing {data_type} skipped - found files already preprocessed')
        return

    logger.info(f'Preprocessing {data_type} starts')

    # Do checks on annotation (sample annot, metadata, patient_regex)
    sample_annotation_df = prep.check_annot(sample_annotation_file, metadata_annotation, prep.in_metadata, patient_regex)

    preprocessed2_file = os.path.join(results_folder, f'preprocessed_{data_type}2.csv')
    if os.path.exists(preprocessed2_file):
        logger.info(f"Reusing previously generated results for {data_type}: {preprocessed2_file}")
        df = pd.read_csv(preprocessed2_file)
    else:
        df = load_sample_data(results_folder, sample_annotation_df, simsi_folder, raw_data_location, run_simsi, run_lfq, debug,
                              data_type)

        preprocess_function = preprocess_fp
        if data_type == 'pp':
            preprocess_function = preprocess_pp

        # ~5 minutes for pp, ~1 hour for fp, of which 1 hour is LFQ
        # returns dataframe in "wide format", i.e. patients as columns
        df = preprocess_function(df, results_folder, picked_fdr, fasta_file, fdr_num_threads, imputation=imputation,
                                 debug=debug, run_lfq=run_lfq)
        df.to_csv(os.path.join(results_folder, preprocessed2_file), index=False)

    sample_annotation_df = sample_annotation.filter_samples_by_metadata(sample_annotation_df, program, entity, histologic_subtype)
    sample_annotation_df.to_csv('/home/mgander/Atlantic/data/sample_annotation_df.csv')

    filtered_sample_annotation_file = os.path.join(results_folder, 'sample_annot_filtered.csv')
    channel_to_sample_id_dict = sample_annotation.get_channel_to_sample_id_dict(sample_annotation_df, filtered_sample_annotation_file, remove_qc_failed=True,
                                                                                remove_replicates=False)
    channel_to_sample_id_dict={a:a.split('Reporter intensity corrected 1 ')[1] for a in df.columns if a[:8]=='Reporter'}

    index_cols = utils.get_index_cols(data_type)
    df_with_ref = prep.rename_columns_with_sample_ids(df, channel_to_sample_id_dict, index_cols=index_cols, remove_ref=False)

    # TODO: should be unneccessary for us? but maybe still a good idea?
    # df_with_ref = prep.remove_ref_empty_batch(df_with_ref, sample_annotation_df)

    df_with_ref.to_csv(os.path.join(results_folder, f'preprocessed_{data_type}_with_ref.csv'), index=False)
    if debug:
        df.to_csv(os.path.join(results_folder, f'debug_preprocessed_{data_type}_before_filter_sample.csv'), index=False)
    df = prep.rename_columns_with_sample_ids(df, channel_to_sample_id_dict, index_cols=index_cols)

    df.to_csv(os.path.join(results_folder, f'preprocessed_{data_type}.csv'), index=False)


def load_sample_data(results_folder: Union[str, Path],
                     sample_annotation_df: pd.DataFrame,
                     simsi_folder: Union[str, Path],
                     raw_data_location: Union[str, Path],
                     run_simsi: bool,
                     run_lfq: bool,
                     debug: bool,
                     data_type: str) -> pd.DataFrame:
    evidence_files = prep.get_evidence_files(sample_annotation_df, raw_data_location, data_type)
    data_loader = TMTLoader(evidence_files)
    if run_simsi:
        data_loader = SimsiTMTLoader(evidence_files, results_folder, simsi_folder, data_type)
    elif run_lfq:
        data_loader = LFQLoader(evidence_files)
    df = prep.load_and_normalize(data_loader, results_folder, data_type=data_type, debug=debug)
    return df


def preprocess_pp(df: pd.DataFrame,
                  results_folder: Union[str, Path],
                  picked_fdr: float,
                  fasta_file: str,
                  fdr_num_threads: int,
                  imputation: bool,
                  debug: bool,
                  run_lfq: bool) -> pd.DataFrame:
    logger.info('Preprocess_pp function')

    # Re-map gene names based on uniprot identifiers in a fasta file. This is necessary 
    # because MaxQuant uses their own uniprot->gene mapping file that cannot be changed.
    df = picked.remap_gene_names(df, fasta_file)

    # create columns to store metadata about the identifications, e.g. imputed, detected in batch, single peptide id
    df = id_meta.create_metadata_columns(df)

    if imputation:
        # Impute missing values within batches
        df.to_csv(os.path.join(results_folder, 'preprocessed_pp_before_imputation.csv'), index=False)
        df = prep.impute_data(df)
        if debug:
            df.to_csv(os.path.join(results_folder, 'debug_preprocessed_pp_after_imputation.csv'), index=False)

    # Filter out contaminants, reverse sequences and non-phospho peptides
    df = prep.filter_data(df, data_type='pp')

    # Aggregate p-peptide intensities across fractions and charge states
    df = prep.sum_peptide_intensities(df, debug, run_lfq)
    if debug:
        df.to_csv(os.path.join(results_folder, 'debug_preprocessed_pp_after_aggregation.csv'), index=False)

    # log10 transform intensities and turn missing values into nans
    df = prep.log_transform_intensities(df)

    # convert to wide format, i.e. each column is a patient with its peptide abundances
    df = prep.convert_long_to_wide_format(df, has_metadata_cols=True)

    # I think solution is to save columns of transfer as separate file and only take to use for reports
    # at least for now
    # test with not running simsi
    if df.columns.str.startswith('Transferred spectra count').any():
        df.loc[:, df.columns.str.startswith('Transferred spectra count')].to_csv(
            os.path.join(results_folder, 'Transfer_metadata.csv'))
        df = df.drop(df.loc[:, df.columns.str.startswith('Transferred spectra count')].columns, axis=1)
    return df


def preprocess_fp(df: pd.DataFrame,
                  results_folder: Union[str, Path],
                  picked_fdr: float,
                  fasta_file: str,
                  fdr_num_threads: int,
                  imputation: bool,
                  debug: bool,
                  run_lfq: bool) -> pd.DataFrame:
    """
    Function for preprocessing full proteome MQ evidence files or simsi results
    :param results_folder: location to save results from picked group fdr + log file
    :param picked_fdr: FDR value for filtering after picked group FDR
    :return: Dataframe with combined preprocessed data from all batches
    """
    logger.info('Preprocess fp function')

    # Apply picked protein group on gene level and filter at 1% FDR
    df = picked.picked_protein_grouping(df, results_folder, picked_fdr, fasta_file, fdr_num_threads)

    # Filter out contaminants and reverse sequences
    df = prep.filter_data(df, data_type='fp')

    # create columns to store metadata about the identifications, e.g. imputed, detected in batch, single peptide id
    df = id_meta.create_metadata_columns(df)

    # Mark proteins detected in the batch but not in the sample
    df = id_meta.mark_detected_in_batch(df)

    # Mark number of peptide identifications per sample
    df = id_meta.mark_num_peptides(df)

    # log10 transform intensities and turn missing values into nans
    df = prep.log_transform_intensities(df)

    return df


if __name__ == '__main__':
    import argparse

    from . import config

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True,
                        help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)

    preprocess_raw(configs["results_folder"], configs["sample_annotation"], configs["simsi"]["simsi_folder"],
                   **configs["preprocessing"], data_types=configs["data_types"])
