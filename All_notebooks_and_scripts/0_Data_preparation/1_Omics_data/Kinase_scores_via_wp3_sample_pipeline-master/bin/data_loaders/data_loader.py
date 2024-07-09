import sys
import re
import logging
from typing import List, Union, Dict
from pathlib import Path

import pandas as pd

from ..utils import get_tmt_channels
from ..preprocess_tools import merge_ms1_columns, median_centering

logger = logging.getLogger(__name__)


class DataLoader:
    def load_data(self, use_cols: List[str]):
        pass

    def median_centering_within_batch(self, all_batches: List[pd.DataFrame]) -> List[pd.DataFrame]:
        dfs = []
        for df in all_batches:
            tmt_channels = get_tmt_channels(df).columns
            df.loc[:, tmt_channels] = median_centering(df.loc[:, tmt_channels])
            dfs.append(df)
        return dfs
    
    def median_centering_ms1(self, all_batches: List[pd.DataFrame]) -> List[pd.DataFrame]:
        """
        Normalizes samples by multiplying each batch MS1 with its own correction factor.
        Only uses peptides detected in >70% of the batches for computing median to prevent
        low-abundant peptides from dragging down the median in runs with deeper coverage
        :param all_batches: list of evidence dataframes
        :return:
        """
        merged_ms1_df = merge_ms1_columns(all_batches)

        # TODO: Find better solution for selecting peptides for median centering. 
        #       Currently, Sarcoma_Batch20-22 have a z-score bias because they have fewer of the "common" peptides than other batches.
        medians = merged_ms1_df[merged_ms1_df.count(axis=1) > 0.7 * len(merged_ms1_df.columns)].median(axis=0).to_dict()
        mean_median = pd.Series(medians.values()).mean()

        dfs = []
        for df in all_batches:
            batch_name = df['Batch'].iloc[0]
            correction_factor = (mean_median / medians[batch_name])
            logger.info(f"Correction factor for {batch_name}: {round(correction_factor, 3)}")

            df['MS1'] = df['Intensity'] * correction_factor
            dfs.append(df)

        return dfs
    
    def impute_ms1_intensity(self, df: pd.DataFrame) -> pd.DataFrame:
        return df
    
    def redistribute_ms1_intensity(self, df: pd.DataFrame) -> pd.DataFrame:
        return df


def extract_cohort_name(evidence_file_path: Union[str, Path]) -> str:
    """Extract batch name including cohort from a file path, e.g. 
    '/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt' => Sarcoma_Batch1
    """
    match = re.search(r'([^/]*)/Batch', evidence_file_path)
    return match.group(1).replace('/', '_')


def extract_batch_name(evidence_file_path: Union[str, Path]) -> str:
    """Extract batch name including cohort from a file path, e.g. 
    '/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt' => Sarcoma_Batch1
    """
    # match = re.search(r'[^/]*/Batch\d+', evidence_file_path)
    match = re.search(r'[^/]*/Batch[^_]+', evidence_file_path)
    return match.group(0).replace('/', '_')


def extract_experiment_name(evidence_file_path: Union[str, Path]) -> str:
    # match = re.search(r'/(Batch\d+_[^/]*)/', evidence_file_path)
    match = re.search(r'/(Batch[\w\d]*_[^/]*)/', evidence_file_path)
    return match.group(1)


def test_batch_names_equals(experiment_to_batch_name_dict: Dict, df: pd.DataFrame):
    # explain
    evidence_files_batch_names = set(experiment_to_batch_name_dict.keys())
    simsi_output_batch_names = set(df['Experiment'].unique().tolist())

    if evidence_files_batch_names != simsi_output_batch_names:
        missing_in_simsi = evidence_files_batch_names - simsi_output_batch_names
        missing_in_evidence_files = simsi_output_batch_names - evidence_files_batch_names
        raise ValueError(f'Batches in list of evidence files and SIMSI evidence do not match:\nMissing in SIMSI: {missing_in_simsi}\nMissing in evidence_files: {missing_in_evidence_files}')
