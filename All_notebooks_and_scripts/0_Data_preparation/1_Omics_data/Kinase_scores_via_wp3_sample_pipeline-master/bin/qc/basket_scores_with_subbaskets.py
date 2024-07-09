import sys
import os

import pandas as pd
import numpy as np

from .. import basket_scoring


"""
!!!! This functionality has been moved to bin/basket_scoring.py !!!!
"""

def main(argv):
    #results_folder = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.12.05_AhS_all_mixed_cohort"
    #results_folder = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.12.14_CJ_batch70"
    results_folder = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.11.24_MT_minimal_test_new_metadata_annot"
    
    metadata_file = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/MasterInform_Metadata_MTB_Portal_221129.xlsx"

    basket_file = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/AS/Pathway Scoring_4th gen/Baskets_4th gen_221214.xlsx"

    all_baskets_df = pd.read_excel(basket_file)
    all_baskets_df['subbasket_level'] = all_baskets_df['SUBBASKET'] + " - " + all_baskets_df['LEVEL']
    all_baskets_df['WEIGHT'] = all_baskets_df['WEIGHT'].fillna(1) # empty cell in WEIGHT column means weight = 1

    z_scores_fp_df = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_z.tsv'), sep='\t')
    z_scores_pp_df = pd.read_csv(os.path.join(results_folder, 'phospho_measures_z.tsv'), sep='\t')

    z_scores_fp_df = z_scores_fp_df.rename(columns=lambda x: x.replace("zscore_", ""))
    z_scores_pp_df = z_scores_pp_df.rename(columns=lambda x: x.replace("zscore_", ""))

    z_scores_pp_df['Modified sequence'] = z_scores_pp_df['Modified sequence'].str.replace(r'([STY])\(Phospho \(STY\)\)', lambda pat: f'p{pat.group(1)}', regex=True)

    protein_phosphorylation_df = pd.read_csv(os.path.join(results_folder, 'protein_results/protein_scores.tsv'), sep='\t')
    kinase_scores_df = pd.read_csv(os.path.join(results_folder, 'kinase_results/kinase_scores.tsv'), sep='\t')

    total_basket_scores = {}
    for basket_name, basket_df in all_baskets_df.groupby('BASKET'):
        print(f"Calculating basket scores for {basket_name}")
        subbasket_scores = {}
        for subbasket_name, subbasket_df in basket_df.groupby('subbasket_level'):
            scoring_rule = subbasket_df['SCORING RULE'].iloc[0].lower()
            print(f"- Calculating subbasket scores for {basket_name} - {subbasket_name} (scoring rule: {scoring_rule})")
            if scoring_rule == 'highest z-score':
                subbasket_score = get_weighted_max(z_scores_fp_df, subbasket_df, 'Gene names', 'GENE NAME')
            elif scoring_rule == 'highest z-score (p-site)':
                subbasket_score = get_weighted_max(z_scores_pp_df, subbasket_df, 'Modified sequence', 'MODIFIED SEQUENCE')
            elif scoring_rule == 'highest protein phosphorylation score (2nd level z-score, fh)':
                subbasket_score = get_weighted_max(protein_phosphorylation_df, subbasket_df, 'Gene names', 'GENE NAME')
            elif scoring_rule == 'highest kinase score (2nd level z-score, fh)':
                subbasket_score = get_weighted_max(kinase_scores_df, subbasket_df, 'PSP Kinases', 'GENE NAME')
            else:
                raise ValueError(f"Unknown scoring rule {scoring_rule}")
            subbasket_scores[f'{basket_name} - {subbasket_name}'] = subbasket_score

        subbasket_scores_df = pd.DataFrame.from_dict(subbasket_scores)

        subbasket_scores_df['{basket_name}_total_basket_score'] = subbasket_scores_df.sum(axis=1)
        total_basket_scores[basket_name] = subbasket_scores_df['{basket_name}_total_basket_score']
        
        # Remove replicate suffix
        subbasket_scores_df['Sample name'] = subbasket_scores_df.index.str.replace(r'-R\d{1}', '', regex=True)

        subtype_df = metadata_df[['Sample name', 'Sarcoma Subtype']]
        subbasket_scores_df = subbasket_scores_df.reset_index().merge(subtype_df, on='Sample name', how='left')
        
        subbasket_output_file = os.path.join(results_folder, f'subbasket_scores_{basket_name}.tsv')
        subbasket_scores_df.to_csv(subbasket_output_file, sep='\t', index=False)
        
        print(f"Written basket results for {basket_name} to: {subbasket_output_file}")
    
    basket_scores_df = pd.DataFrame.from_dict(total_basket_scores)
    return basket_scores_df


def get_weighted_max(z_score_df: pd.DataFrame, subbasket_df: pd.DataFrame, z_score_index_col: str, subbasket_index_col: str):
    z_scores = get_weighted_z_scores(z_score_df, subbasket_df, z_score_index_col, subbasket_index_col)
    
    # take the maximum score per column (=sample)
    return z_scores.max()


def get_weighted_z_scores(z_score_df: pd.DataFrame, subbasket_df: pd.DataFrame, z_score_index_col: str, subbasket_index_col: str):
    """
    z_score_df is a pandas dataframe with z_scores as columns and genes/modified sequences as rows
    subbasket_df is a pandas dataframe with genes with weights as rows
    z_score_index_col is the column name with the identifier in z_score_df, e.g. "Gene names" or "Modified sequence"
    subbasket_index_col is the column name with the identifier in subbasket_df, e.g. "GENE NAME" or "Modified sequence"
    """
    # only keep samples and z_score_index_col columns
    z_scores = z_score_df.set_index(z_score_index_col)
    z_scores = z_scores.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+') 
    z_scores = z_scores.astype(float) # make sure all z-scores are floats and not objects
    z_scores = z_scores.reset_index()
    
    # deal with protein groups; in the z-score dataframe protein groups exists, 
    # e.g. AKT1;AKT2;AKT3. We temporarily split them in separate rows to merge
    # with subbasket_df, where each row only contains a single protein. After 
    # merging we combine those rows into a single row again.
    z_scores, z_score_index_col_exploded = explode_on_separated_string(z_scores, z_score_index_col)
    
    # merge z-score dataframe with basket genes
    z_scores = z_scores.merge(subbasket_df[[subbasket_index_col, 'WEIGHT']], left_on=z_score_index_col_exploded, right_on=subbasket_index_col)
    if set(subbasket_df[subbasket_index_col]) != set(z_scores[z_score_index_col_exploded]):
        print(f"WARNING: could not find all identifiers in the z_score data: {set(subbasket_df[subbasket_index_col])-set(z_scores[z_score_index_col_exploded])}")
    
    # merge proteins from the same protein group into a single row
    z_scores = z_scores.groupby(z_score_index_col).agg('first').reset_index()
    
    # multiply the z-score by the associated weight (usually -1 or +1)
    weights = z_scores['WEIGHT']
    z_scores = z_scores.drop(columns=[z_score_index_col_exploded, z_score_index_col, subbasket_index_col, 'WEIGHT'])
    z_scores = z_scores.multiply(weights, axis=0)
    
    return z_scores


def explode_on_separated_string(df: pd.DataFrame, index_col: str, sep: str=';'):
    index_col_exploded = f'{index_col}_exploded'
    df[index_col_exploded] = df[index_col].str.split(sep)
    return df.explode(index_col_exploded), index_col_exploded


if __name__ == "__main__":
    main(sys.argv[1:])
