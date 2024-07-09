import sys
from typing import Dict
import logging

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from .. import sample_annotation
from .. import identification_metadata as id_meta
from .. import utils
from .. import basket_scoring

pd.set_option('display.max_rows', None)

# hacky way to get the package logger instead of just __main__ when running as a module
logger = logging.getLogger(__package__ + "." + __file__)


def plot_batch_overlaps(results_folder: str, sample_annotation_file: str, baskets_file: str, data_type: str, group_by_batch: bool=False, plot_stable_detections: bool=False, plot_proteins_of_interest: bool=False):
    sample_annotation_df = sample_annotation.load_sample_annotation(sample_annotation_file)
    channel_to_sample_id_dict = sample_annotation.get_channel_to_sample_id_dict(sample_annotation_df)
    
    # preprocessed_{fp,pp}2.csv still contain batch numbers in the column names, as well as protein information such as number of detected peptides
    preprocessed2_file = f'{results_folder}/preprocessed_{data_type}2.csv'
    
    index_col = 'Protein IDs'
    y_label = 'proteins'
    if data_type == 'pp':
        index_col = 'Modified sequence'
        y_label = 'phosphopeptides'
    
    keep_cols = [index_col]
    units = 'samples'
    if group_by_batch:
        units = 'batches'
    
    figure_suffix = ""
    
    df = pd.read_csv(preprocessed2_file, 
                     index_col=index_col, 
                     usecols=lambda x: x.startswith('Reporter intensity corrected') or x in keep_cols)
    df.columns = utils.remove_cohort_names_from_channel_names(df.columns)
    df = df.filter(items=channel_to_sample_id_dict.keys(), axis=1)
    
    plot_overlap_histogram(df, group_by_batch)
    
    labels = ['all proteins']
    if plot_stable_detections:
        if data_type == "pp":
            raise ValueError("Currently, no definition of stability available for phospho data")
        figure_suffix = "_stable"
        
        # this looks at number of peptides across all samples
        annot_df = pd.read_csv(preprocessed2_file, 
                         index_col=index_col, 
                         usecols=lambda x: x.startswith(id_meta.METADATA_COLUMN_PREFIX) or x in keep_cols)
        annot_df.columns = utils.remove_cohort_names_from_channel_names(annot_df.columns)
        stable_df = id_meta.filter_by_min_peptides(df.rename(columns=lambda x: x.replace('Reporter intensity corrected ', '')), annot_df, 2)
        
        plot_overlap_histogram(stable_df, group_by_batch)
        labels += ['>= 2 peptides']
    
    if plot_proteins_of_interest:
        figure_suffix += "_proteins_of_interest"
        all_baskets_df = basket_scoring.read_baskets_file_4th_gen(baskets_file)
        basket_index_col = 'GENE NAME'
        scoring_rule = 'highest z-score'
        if data_type == "pp":
            scoring_rule = 'highest z-score (p-site)'
            basket_index_col = 'MODIFIED SEQUENCE'
            
        proteins_of_interest_df = all_baskets_df[all_baskets_df['SCORING RULE'] == scoring_rule][[basket_index_col, 'WEIGHT']]
        
        df = df.reset_index()
        df, df_index_col_exploded = utils.explode_on_separated_string(df, index_col)
        df = df.merge(proteins_of_interest_df, left_on=df_index_col_exploded, right_on=basket_index_col, how='right')
        if set(proteins_of_interest_df[basket_index_col]) != set(df[df_index_col_exploded]):
            logger.warning(f"Could not find all identifiers: {set(proteins_of_interest_df[basket_index_col])-set(df[df_index_col_exploded])}")
        
        df = df.drop(columns=[basket_index_col, df_index_col_exploded, 'WEIGHT']).set_index(index_col)
        
        plot_overlap_histogram(df, group_by_batch)
        labels += ['proteins of interest']
    
    
    plt.ylabel(f'#{y_label}')
    plt.xlabel(f'detected in at least n {units}')
    plt.legend(labels)
    plt.tight_layout()
    #plt.show()

    output_file = f'{results_folder}/{y_label}_{units}_overlaps_{data_type}{figure_suffix}.jpg'
    plt.savefig(output_file, bbox_inches='tight', dpi=400)
    
    output_file = f'{results_folder}/{y_label}_{units}_overlaps_{data_type}{figure_suffix}.eps'
    plt.savefig(output_file, bbox_inches='tight')


def plot_overlap_histogram(df: pd.DataFrame, group_by_batch: bool):    
    # transpose dataframe such that samples are rows
    intensity_df = df.T
    intensity_df = intensity_df.reset_index()
    
    if group_by_batch:
        intensity_df['Batch'] = intensity_df['index'].apply(lambda x: x.split()[-1])
        proteins_by_batch_df = intensity_df.groupby('Batch').agg('sum')
    else:
        proteins_by_batch_df = intensity_df
    
    num_batches = len(proteins_by_batch_df.index)

    num_batches_per_protein_df = proteins_by_batch_df.replace(0, np.nan).count(axis=0)

    num_batches_per_protein_df.hist(cumulative=-1, bins=range(1, num_batches+1))
    
    # round up to the nearest multiple of 50
    plt.xlim([0, (int(num_batches / 50) + 1) * 50])


if __name__ == '__main__':
    """
    python -m bin.qc.get_protein_batch_overlaps
    """
    import argparse
    
    from .. import config
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default='/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.12.29_AhS_batch73_74_imp/configs.json',
                        help="Absolute path to configuration file.")
    parser.add_argument("-d", "--data_type", default='fp',
                        help="fp for full proteome and pp for phosphoproteome.")
    parser.add_argument("-g", "--group_by_batch", default=False, action='store_true',
                        help="Group samples by batch and make plots with number of batches on the x-axis.")
    parser.add_argument("-s", "--plot_stable_detections", default=False, action='store_true',
                        help="Only count proteins that are identified by at least 2 peptides.")
    parser.add_argument("-i", "--plot_proteins_of_interest", default=False, action='store_true',
                        help="Only plot proteins of interest based on the basket annotations.")

    
    args = parser.parse_args(sys.argv[1:])

    configs = config.load(args.config)
    
    configs["clinic_proc"]["prot_baskets"] = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/AS/Pathway Scoring_4th gen/Baskets_4th gen_230103.xlsx"
    
    plot_batch_overlaps(configs["results_folder"], configs["sample_annotation"], configs["clinic_proc"]["prot_baskets"], args.data_type, args.group_by_batch, args.plot_stable_detections, args.plot_proteins_of_interest)

