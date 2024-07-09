import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Union

import pandas as pd

import bin.config as config
import bin.TOPAS_scoring_functions as scoring

logger = logging.getLogger(__name__)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def read_kinase_scoring(results_folder: Union[str, Path]):
    if os.path.exists(os.path.join(results_folder, 'kinase_results', 'kinase_scores.tsv')):
        try:
            kinase_scores = os.path.join(results_folder, 'kinase_results', 'kinase_scores.tsv')
            kinase_scores_df = pd.read_csv(kinase_scores, index_col=['PSP Kinases', 'No. of total targets'], sep='\t')
        except PermissionError:
            raise PermissionError(
                f'Cannot open kinase scores file, check if you have it open in Excel. {kinase_scores}')
    return kinase_scores_df


def explode_psites_by_kinases(patient_df):
    """Filter for p-sites with a kinase annotation and explode by the kinase annotation.

    Args:
        patient_df: dataframe with a 'PSP Kinases' column that can contain multiple kinases.

    Returns:
        Exploded dataframe with a single kinase in the PSP Kinases column
    """
    psites_with_kinases = patient_df['PSP Kinases'].notna()
    patient_df.loc[psites_with_kinases, 'PSP Kinases'] = patient_df.loc[psites_with_kinases, 'PSP Kinases'].str.split(";")
    patient_df = patient_df.loc[psites_with_kinases, :].explode('PSP Kinases')
    return patient_df


def generate_kinase_scores():
    pass


# TODO: keep psite column from psite annotation in kinase peptide dataframe


def kinase_scoring(results_folder, preprocessed_df):
    logger.info('Running kinase scoring module')

    if os.path.exists(os.path.join(results_folder, 'kinase_results', 'kinase_scores.tsv')):
        logger.info(f'Kinase scoring skipped - found files already processed')
        return

    if not os.path.exists(os.path.join(results_folder, 'kinase_results')):
        os.makedirs(os.path.join(results_folder, 'kinase_results'))

    logger.info('  Calculate p-site weights')
    kinase_df = scoring.calculate_psite_weights(preprocessed_df)

    kinase_df = explode_psites_by_kinases(kinase_df)
    kinase_df.to_csv(os.path.join(results_folder, 'kinase_results', 'kinase_annotated_patients_with_weights.tsv'), sep='\t', index=False)

    logger.info('  Calculate modified sequence weights')
    kinase_summed_weights = scoring.calculate_modified_sequence_weights(kinase_df, 'PSP Kinases')

    # TODO: Export here; table with uncapped weights and uncapped zscores
    kinase_capped_values = scoring.cap_zscores_and_weights(kinase_summed_weights)

    logger.info('  Calculate weighted z-scores')
    kinase_scored_peptides = scoring.calculate_weighted_z_scores(kinase_capped_values)
    kinase_scored_peptides.to_csv(os.path.join(results_folder, 'kinase_results', 'scored_peptides.tsv'), sep='\t',
                                  index=False)

    logger.info('  Calculate kinase scores')
    kinase_first_level_scores = scoring.sum_weighted_z_scores(kinase_scored_peptides, by='PSP Kinases')

    # scoring.plot_histograms_to_check_normality(kinase_first_level_scores)

    logger.info('  2nd level z-scoring, adding target space and writing results')
    kinase_scores = scoring.second_level_z_scoring(kinase_first_level_scores, 'PSP Kinases')
    kinase_spaces = scoring.get_target_space(annotated_peptides_df=kinase_df, scored_peptides_df=kinase_scored_peptides, grouping_by='PSP Kinases')
    kinase_scores = pd.merge(left=kinase_spaces, right=kinase_scores, on='PSP Kinases', how='left').sort_values(
        by='PSP Kinases')
    kinase_scores = kinase_scores.set_index(['PSP Kinases', 'No. of total targets'])
    kinase_scores.to_csv(os.path.join(results_folder, 'kinase_results', 'kinase_scores.tsv'), sep='\t')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default='/home/fhamood/PycharmProjects/WP3_Pipeline/wp3_sample_pipeline/config_patients.json',
                        help="Absolute path to configuration file.")
    args = parser.parse_args(sys.argv[1:])
    configs = config.load(args.config)

    preprocessed_path = os.path.join(configs["results_folder"], 'topas_score_preprocessed.tsv')
    if os.path.exists(preprocessed_path):
        print('Found preprocessed file')
        preprocessed_df = pd.read_csv(preprocessed_path, sep = '\t')
    else:
        preprocessed_df = scoring.topas_score_preprocess(configs["results_folder"], configs["patient_regex"])

    kinase_scoring(configs["results_folder"], preprocessed_df)
