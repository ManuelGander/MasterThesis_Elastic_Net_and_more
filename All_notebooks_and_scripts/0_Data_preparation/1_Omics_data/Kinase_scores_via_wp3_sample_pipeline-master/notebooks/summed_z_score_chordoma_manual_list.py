# manual cureated signatures of chordoma across all patients
path_to_signauture_file = '/home/amir/Desktop/Annika_files/IFN TUPAC score_CHDM_final.xlsx'

import json
from pathlib import  Path
import pandas as pd


result_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.08.01_CJ_paper_pdx_chordomacl'
meta_path= '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Metadata_Papercohort_CHDMmodels_230801.xlsx'
z_scores_file_path = 'full_proteome_measures_z.tsv'
z_scores_path = Path(result_folder)/Path(z_scores_file_path)


z_scores_df = pd.read_csv(z_scores_path,sep='\t').set_index('Gene names')

signature = pd.read_excel(path_to_signauture_file)
signatures_proteins = signature.iloc[:,0].tolist()

alpha_signatures_df = z_scores_df[z_scores_df.index.isin(signatures_proteins)]
sum_z_scores = pd.DataFrame(alpha_signatures_df.sum(axis=0,skipna=True)).T
sum_z_scores.index = ['sum']
final = pd.concat([alpha_signatures_df,sum_z_scores],axis=0)

meta_df = pd.read_excel(meta_path)
z_scores_df = final.T

# adding LOO z_scores
df = pd.DataFrame(z_scores_df['sum'].copy())
df.columns = ['sum']
new_z_scors=[]
for i in df.index:
    loo_df = df.drop(i, axis=0)
    median = float(loo_df.median())
    stdev = float(loo_df.std())
    new_z_scors.append((df.loc[i,'sum'] - median) / stdev)
z_scores_df['z_scores_sum_LOO'] = new_z_scors


z_scores_df["Sample name"] = z_scores_df.index
z_scores_df["Sample name"] = z_scores_df["Sample name"].str.replace('zscore_','')

z_scores_df = z_scores_df.merge(meta_df,on='Sample name')
z_scores_df.index = z_scores_df["Sample name"]


z_scores_df.to_excel('/home/amir/Desktop/Annika_files/TUPAC_for_new_signatures.xlsx')