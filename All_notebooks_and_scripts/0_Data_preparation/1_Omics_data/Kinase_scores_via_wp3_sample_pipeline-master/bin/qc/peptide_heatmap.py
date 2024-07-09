import sys
import os
import json
import logging
from typing import Dict, List

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.cm as colormap
from matplotlib.colors import ListedColormap
from scipy.stats import zscore
import picked_group_fdr.digest as digest

import bin.config as config
from bin.preprocess_tools import rename_columns_with_sample_ids
from bin.sample_annotation import get_channel_to_sample_id_dict

logger = logging.getLogger(__name__)

pd.set_option('display.max_rows', None)


def plot_heatmaps(df, gene_name, batches_of_interest, channel_to_sample_id_dict, patients_of_interest, fasta_file,
                  protein_group_gene_names, output_folder):
    df_filtered = df[(df['Gene names'].apply(lambda x: gene_name in x.split(';') if str(x) != 'nan' else False)) & (df['Modifications'] == 'Unmodified')]
    #logger.info(df_filtered[['Modified sequence', 'Charge']].value_counts())
    
    tmt_cols = [f'Reporter intensity corrected {i}' for i in range(1,9)]
    index_cols = ['Modified sequence', 'Charge']
    meta_cols = ['PEP', 'MS1']
    summed_df = df_filtered.groupby(['Batch'] + index_cols).agg('sum', numeric_only=True).reset_index()
    indexed_df = summed_df.set_index(index_cols)[['Batch'] + tmt_cols + meta_cols]
    
    wide_df = convert_long_to_wide_format(indexed_df)
    wide_df = wide_df.reset_index()
    wide_df = add_sequence_positions(gene_name, wide_df, fasta_file, protein_group_gene_names)
    best_covered_isoform = list(wide_df.filter(regex='Sequence position').count().sort_values(ascending=False).index)
    wide_df = wide_df.sort_values(by=best_covered_isoform)
    
    plot_peptide_heatmap(gene_name, wide_df, batches_of_interest, channel_to_sample_id_dict, patients_of_interest, index_cols, output_folder)
    plot_PEP_heatmap(gene_name, wide_df, batches_of_interest, output_folder)
    plot_MS1_heatmap(gene_name, wide_df, batches_of_interest, output_folder)
    plot_isoform_heatmap(gene_name, wide_df, output_folder)


def add_sequence_positions(gene_name, wide_df, fasta_file, protein_group_gene_names):
    sequences = list()
    for (protein, _, fasta_gene_name, _), (_, sequence) in zip(digest.readFastaProteins(fasta_file, db='target'), digest.readFastaMaxQuant(fasta_file, db='target')):
        if fasta_gene_name == gene_name or fasta_gene_name in protein_group_gene_names:
            sequences.append((protein, fasta_gene_name, sequence))
            #logger.info(f'>{protein} GN={fasta_gene_name}\n{sequence}')
    for protein, fasta_gene_name, sequence in sequences:
        wide_df[f"Sequence position {fasta_gene_name} {protein}"] = wide_df['Modified sequence'].apply(lambda x: sequence.index(x[1:-1]) if x[1:-1] in sequence else np.nan)
    return wide_df


def plot_isoform_heatmap(gene_name, wide_df, output_folder):
    isoform_df = wide_df.set_index(['Modified sequence', 'Charge']).filter(regex='Sequence position').rename(columns=lambda x: x.replace('Sequence position ', ''))
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(isoform_df.T, xticklabels=True, yticklabels=True, vmin=0, vmax=isoform_df.max().max(), cmap="Greens", cbar_kws={'label': 'Sequence position start'})
    plt.tight_layout()
    plt.savefig(f'{output_folder}/{gene_name}_isoforms.pdf')


def plot_MS1_heatmap(gene_name, wide_df, batches_of_interest, output_folder):
    intensity_df = wide_df.set_index(['Modified sequence', 'Charge'])[[f'MS1 {b}' for b in sorted(list(batches_of_interest))]].apply(np.log10)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(intensity_df.T, xticklabels=True, yticklabels=True, vmin=6, vmax=8, cmap="coolwarm", cbar_kws={'label': 'log10 MS1 intensity'})
    plt.tight_layout()
    plt.savefig(f'{output_folder}/{gene_name}_ms1_intensity.pdf')


def plot_PEP_heatmap(gene_name, wide_df, batches_of_interest, output_folder):
    pep_df = wide_df.set_index(['Modified sequence', 'Charge'])[[f'PEP {b}' for b in sorted(list(batches_of_interest))]].apply(np.log10)
    pep_df = pep_df.replace(-1*np.inf, 0.01)
    
    viridis = colormap.get_cmap('coolwarm', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    darkGrey = np.array([64/256, 64/256, 64/256, 1])
    newcolors[255, :] = darkGrey
    viridis_ext = ListedColormap(newcolors)
    
    plt.figure(figsize=(10,8))
    sns.heatmap(pep_df.T, xticklabels=True, yticklabels=True, vmin=-5, vmax=0, cmap=viridis_ext, cbar_kws={'label': 'log10 error probability'})
    plt.tight_layout()
    plt.savefig(f'{output_folder}/{gene_name}_id_error.pdf')


def plot_peptide_heatmap(gene_name, wide_df, batches_of_interest, channel_to_sample_id_dict, patients_of_interest, index_cols, output_folder):
    channel_to_sample_id_dict_filtered = {k: f'{v}-{k.split()[-1]}' for k, v in channel_to_sample_id_dict.items() if k.split()[-1] in batches_of_interest}
    peptide_df = rename_columns_with_sample_ids(wide_df, channel_to_sample_id_dict_filtered, index_cols=index_cols)
    peptide_df = peptide_df.set_index(index_cols)
    peptide_df = peptide_df.rename(columns=lambda x: f"{x}-POI" if '-'.join(x.split('-')[:-1]) in patients_of_interest else x)
    zscore_df = peptide_df.replace(0, np.nan).apply(np.log10)
    
    zscore_df = zscore_df.sub(zscore_df.mean(axis=1), axis=0)
    zscore_df = zscore_df.div(zscore_df.std(axis=1), axis=0)
    
    plt.figure(figsize=(18, 16))
    sns.heatmap(zscore_df.T, cmap="coolwarm", xticklabels=True, yticklabels=True)
    plt.tight_layout()
    plt.savefig(f'{output_folder}/{gene_name}_z_scores.pdf')


def convert_long_to_wide_format(indexed_df):
    all_batches = list()
    for batch_name, batch_df in indexed_df.groupby('Batch'):
        batch_df = batch_df.drop(columns='Batch')
        batch_df = batch_df.rename(columns=lambda x: f'{x} {batch_name}')
        all_batches.append(batch_df)
    wide_df = pd.DataFrame().join(all_batches, how="outer")
    return wide_df


def plot_peptide_heatmaps(results_folder, sample_annotation):
    # gene_name = 'IGFBP5'
    # gene_name = 'PTEN'
    # protein_group_gene_names = []
    gene_name = 'CDKN2A'
    protein_group_gene_names = ['CDKN2B']
    # gene_name = 'NTRK2'
    # protein_group_gene_names = []

    # patients_of_interest = ['H021-1KJALE-M1-E2','H021-JLP1D4-T2','H021-RCJFWM-T1','H021-RF6GGV-M1','H021-RF6GGV-T1','H021-S5LP5E-M2-E2','H021-W834UD-T2-E2','H021-XQGYFN-T2-E2']
    patients_of_interest = []

    df = pd.read_csv(os.path.join(results_folder, 'debug_preprocessed_fp_after_ms1_correction.csv'))

    fasta_file = '/media/kusterlab/internal_projects/active/TOPAS/Databases/uniprot_proteome_up000005640_03112020.fasta'

    sample_annotation_df = pd.read_csv(sample_annotation)
    channel_to_sample_id_dict = get_channel_to_sample_id_dict(sample_annotation_df, remove_qc_failed=True)

    sample_id_to_channel_dict = {v: k for k, v in channel_to_sample_id_dict.items()}
    batches_of_interest = {sample_id_to_channel_dict[p].split()[-1] for p in patients_of_interest}
    
    all_batches = {sample_id.split()[-1] for sample_id in channel_to_sample_id_dict.keys()}
    if len(all_batches) < 5:
        batches_of_interest |= all_batches
    else:
        batches_of_interest |= {'Sarcoma_Batch18', 'Sarcoma_Batch23', 'Sarcoma_Batch29', 'Sarcoma_Batch32', 'Sarcoma_Batch50',
                                'Sarcoma_Batch53', 'Sarcoma_Batch55', 'Sarcoma_Batch57'}
    logger.info(f'Batches of interest: {batches_of_interest}')

    plot_heatmaps(df, gene_name, batches_of_interest, channel_to_sample_id_dict, patients_of_interest, fasta_file, protein_group_gene_names,
                  results_folder)


if __name__ == "__main__":
    config_file = sys.argv[1]
    configs = config.load(config_file)

    plot_peptide_heatmaps(configs["results_folder"], configs['sample_annotation'])
