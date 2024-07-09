import re
import sys
from pathlib import Path
import logging
import requests
import time
import json
from typing import List

import pandas as pd
import numpy as np


def get_phospho_score_and_related_peptides_z_scores(peptides_z_scores_df: pd.DataFrame,
                                                    ps_scores_df: pd.DataFrame,
                                                    protein_of_interest: str,
                                                    patients_list: list) -> pd.DataFrame:
    """ 
    Gets the Phosphorylation scores and peptides z_scores given a protein and list of patients
    for figure 2C of the MTB paper
    :peptides_z_scores_df:  phospho_measures_z.tsv
    :ps_scores_df: protein_scores.tsv
    :protein_of_interest: like EGFR
    :patients_list: a list of patients
    """
    ps_scores_df = ps_scores[ps_scores['Gene names'] == protein_of_interest]
    peptides_z_scores_df = phospho_measure[phospho_measure['Gene names'] == protein_of_interest]
    peptides_z_scores_df = peptides_z_scores_df.set_index('Modified sequence')
    ps_scores_df = ps_scores_df.set_index('Gene names')

    peptides_z_scores_df = peptides_z_scores_df.filter(regex='zscore_')
    peptides_z_scores_df.columns = peptides_z_scores_df.columns.str.replace('zscore_', '')
    peptides_z_scores_df = peptides_z_scores_df[patients_list]

    ps_scores_df = ps_scores_df[patients_list]
    if all(ps_scores_df.columns == peptides_z_scores_df.columns):
        final_df = pd.concat([peptides_z_scores_df, ps_scores_df])
    return final_df


def get_kinase_score_and_related_peptides_scores(pep_kinase_Scores_df: pd.DataFrame,
                                                 kinase_scores_df: pd.DataFrame,
                                                 kinase_of_interest: str,
                                                 patients_list: list) -> pd.DataFrame:
    """ 
    Gets the kinase scors and peptides scores from the kinase for a kinase and list of patients
    for figure 2D of MTB paper
    kinase_scores.tsv  scored_peptides.tsv
    :pep_kinase_Scores_df:  scored_peptides.tsv
    :kinase_scores_df: kinase_scores.tsv
    :kinase_of_interest: like EGFR
    :patients_list: a list of patients
    """
    pep_kinase_Scores_df = pep_kinase_Scores_df[pep_kinase_Scores_df['PSP Kinases'] == kinase_of_interest]
    kinase_scores_df = kinase_scores_df[kinase_scores_df['PSP Kinases'] == kinase_of_interest]
    pep_kinase_Scores_df = pep_kinase_Scores_df.set_index('Modified sequence')
    kinase_scores_df = kinase_scores_df.set_index('PSP Kinases')

    pep_kinase_Scores_df = pep_kinase_Scores_df.filter(regex='weighted_')
    pep_kinase_Scores_df.columns = pep_kinase_Scores_df.columns.str.replace('weighted_', '')
    pep_kinase_Scores_df = pep_kinase_Scores_df[patients_list]

    kinase_scores_df = kinase_scores_df[patients_list]
    if all(kinase_scores_df.columns == pep_kinase_Scores_df.columns):
        final_df = pd.concat([pep_kinase_Scores_df, kinase_scores_df])
    return final_df


def get_pep_number_from_protein_name(intensity_df: pd.DataFrame, protein_name: str) -> pd.DataFrame:
    """ 
    gets the number of the pepetides from the protein name for all the patients
    :intensity_df: A pandas dataframe of the intensities with one 'Identification metadata' column per  patient; from preprocessing step 
    :protein_name: the name of the protein
    :USAGE :
        get_pep_number_from_protein_name(fp_df,'EGFR')

    """
    try:
        protein_df = pd.DataFrame(intensity_df.loc[protein_name, :])
        premeta_df = protein_df.filter(regex='^Identification metadata ', axis=0)
        premeta_df[protein_name] = pd.to_numeric(premeta_df[protein_name].str.replace('num_peptides=|;', ''))
        premeta_df['Sample name'] = premeta_df.index.str.replace('Identification metadata ', '')
        premeta_df.columns = ['num_pep', 'Sample name']
        return premeta_df
    except:
        pass


def get_peptides_for_protein_of_interest_from_evidence_file(protein_interest: str,
                                                            evidence: pd.DataFrame,
                                                            sample_annotation_df: pd.DataFrame,
                                                            TMT_no=11) -> pd.DataFrame:
    """
    Returns the peptides related to one protein from the evidence file as a DataFrame

    :protein_interest: the protein name to extract the peptides ; only one protein 
    :evidence: the evidence file from MaxQuant
    :sample_annotation_df: meta_data file with the columns 'Batch Name'   'TMT channel' and 'Cohort'
    :TMT_no: number of channels/plexes 

    :USAGE: 
            get_peptides_for_protein_of_interest_from_evidence_file('EGFR',pp_df,sample_annotation_df)

    """

    sample_annotation_df = sample_annotation_df[sample_annotation_df['QC'] != 'failed']
    sub_set_df = evidence[evidence['Gene names'] == protein_interest]
    list_cols = [*['Modified sequence', 'Batch'], *[f'Reporter intensity corrected {i + 1}' for i in range(TMT_no)]]
    sub_set_df = sub_set_df[list_cols]
    long_df = pd.melt(sub_set_df,
                      id_vars=['Modified sequence', 'Batch'],
                      value_vars=[f'Reporter intensity corrected {i + 1}' for i in range(TMT_no)])
    long_df["sample"] = long_df['variable'] + '_' + long_df['Batch']
    sample_annotation_df['sample'] = 'Reporter intensity corrected ' + sample_annotation_df["TMT Channel"].astype(str) + '_' + \
                                     sample_annotation_df["Cohort"] + '_Batch' + sample_annotation_df["Batch Name"].astype(str)

    final_df = sample_annotation_df.merge(long_df, on='sample')
    return final_df


def init_file_logger(results_folder, log_file_name):
    module_name = ".".join(__name__.split(".")[:-1])
    file_logger = logging.FileHandler(results_folder / Path(log_file_name))
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    file_logger.setFormatter(formatter)
    logging.getLogger(module_name).addHandler(file_logger)


def send_slack_message(message):
    slack_data = {
        "username": "WP3 sample pipeline",
        "icon_emoji": ":gem:",
        "channel": "#wp3_pipeline",
        "attachments": [{"fields": [{"title": 'New Incoming Message', "value": message}]}]
    }
    webhook_url = "https://hooks.slack.com/services/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"
    response = requests.post(
        webhook_url, data=json.dumps(slack_data),
        headers={'Content-Type': 'application/json', 'Content-Length': str(sys.getsizeof(slack_data))}
    )
    if response.status_code != 200:
        raise ValueError(
            'Request to slack returned an error %s, the response is:\n%s'
            % (response.status_code, response.text)
        )


def get_index_cols(data_type: str) -> List[str]:
    index_cols = ['Gene names']
    if data_type == 'pp':
        index_cols = ['Gene names', 'Modified sequence', 'Proteins']
    return index_cols


def short_phospho_notation(modified_sequence_column: pd.Series) -> pd.Series:
    return modified_sequence_column.str.replace(r'([STY])\(Phospho \(STY\)\)', lambda pat: f'p{pat.group(1)}', regex=True)


def split_str_list(samples_for_report):  # add typehint
    samples_for_report = re.split(r'[;|,]', samples_for_report)
    samples_for_report = [sample.strip() for sample in samples_for_report]
    return samples_for_report


def whitespace_remover(df):
    for col in df.columns:
        if df[col].dtype == 'object':
            # applying strip function on column
            df[col] = df[col].astype(str).map(str.strip)
            df[col] = df[col].replace(r'^nan$', np.nan, regex=True)
        else:
            pass
    return df


def remove_cohort_names_from_channel_names(index: pd.Index):
    return index.str.replace(r'(\d{1,2}) .*(Batch\d{1,3})', r'\1 \2', regex=True)


def get_tmt_channels(df: pd.DataFrame) -> pd.DataFrame:
    return df.filter(regex=r'^Reporter intensity corrected \d{1,2}')


def keep_only_sample_columns(df: pd.DataFrame, patient_regex: str) -> pd.DataFrame:
    # Filter for patient expression columns, sample IDs have the format 
    # - "I007-031-108742" (INFORM)
    # - "H021-7AAYWW-T1" (MASTER/Chordoma)
    # - "P024567" (Glioma)
    # - "Reporter intensity corrected" (QC channels)
    # TODO: move regex to config file
    return df.filter(regex=rf'({patient_regex})|(Reporter intensity corrected)|(^P\d{6}$)')


def explode_on_separated_string(df: pd.DataFrame, index_col: str, sep: str = ';'):
    index_col_exploded = f'{index_col}_exploded'
    df[index_col_exploded] = df[index_col].str.split(sep)
    return df.explode(index_col_exploded), index_col_exploded
