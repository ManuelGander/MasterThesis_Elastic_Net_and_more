#https://pubmed.ncbi.nlm.nih.gov/30559422/
paper_sigs = '/home/amir/Desktop/Annika_files/signatures/paper.txt'  # https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-018-0302-5/MediaObjects/41591_2018_302_MOESM2_ESM.xlsx


import json
from pathlib import  Path
import pandas as pd


result_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.08.01_CJ_paper_pdx_chordomacl'
meta_path= '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Metadata_Papercohort_CHDMmodels_230801.xlsx'
z_scores_file_path = 'full_proteome_measures_z.tsv'
z_scores_path = Path(result_folder)/Path(z_scores_file_path)



def unnest_proteingroups(df:pd.DataFrame) -> pd.DataFrame:
    """
    Unnest the protein_groups A;B as two separate rows with the same values
    the protein groups are the index of the the pandas dataframe df
    """
    temp_df = df
    temp_df['index'] = temp_df.index.str.split(';')
    temp_df = temp_df.explode('index')
    temp_df = temp_df.set_index('index')
    return temp_df



signatures_proteins = pd.read_csv(paper_sigs)
signatures_proteins = signatures_proteins.iloc[:,0].tolist()

z_scores_df = pd.read_csv(z_scores_path,sep='\t').set_index('Gene names')
z_scores_df = unnest_proteingroups(z_scores_df)
alpha_signatures_df = z_scores_df[z_scores_df.index.isin(signatures_proteins)]
sum_z_scores = pd.DataFrame(alpha_signatures_df.sum(axis=0,skipna=True)).T
sum_z_scores.index = ['sum']
final = pd.concat([alpha_signatures_df,sum_z_scores],axis=0)

meta_df = pd.read_excel(meta_path)
z_scores_df = final.T
z_scores_df["Sample name"] = z_scores_df.index
z_scores_df["Sample name"] = z_scores_df["Sample name"].str.replace('zscore_','')
z_scores_df = z_scores_df.merge(meta_df,on='Sample name')
z_scores_df.index = z_scores_df["Sample name"]
z_scores_df.to_excel('~/Desktop/nature_paper_signatures.xlsx')