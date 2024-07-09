import json
from pathlib import  Path
import pandas as pd

# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.html
alpha_signatures_path = '/home/amir/Desktop/Annika_files/HALLMARK_INTERFERON_ALPHA_RESPONSE.v2023.1.Hs.json'

# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_INTERFERON_GAMMA_RESPONSE.html
#alpha_signatures_path = '/home/amir/Desktop/Annika_files/HALLMARK_INTERFERON_GAMMA_RESPONSE.v2023.1.Hs.json'


result_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.08.01_CJ_paper_pdx_chordomacl'
meta_path= '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Metadata_Papercohort_CHDMmodels_230801.xlsx'
z_scores_file_path = 'full_proteome_measures_z.tsv'
z_scores_path = Path(result_folder)/Path(z_scores_file_path)


def json_reader(config_path):
    with open(config_path, 'r') as f:
        data = json.load(f)
    return data



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



signatures_proteins = json_reader(alpha_signatures_path)['HALLMARK_INTERFERON_ALPHA_RESPONSE']['geneSymbols']
#signatures_proteins = json_reader(alpha_signatures_path)['HALLMARK_INTERFERON_GAMMA_RESPONSE']['geneSymbols']

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
z_scores_df.to_excel('~/Desktop/alpha_signatures_paper+paper_moded.xlsx')