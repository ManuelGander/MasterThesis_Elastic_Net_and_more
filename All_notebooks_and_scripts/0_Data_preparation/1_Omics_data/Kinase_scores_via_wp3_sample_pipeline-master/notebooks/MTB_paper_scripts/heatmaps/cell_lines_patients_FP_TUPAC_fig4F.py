
# heatmap plot
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
from pathlib import Path

report_dir = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.06.22_AhS_PAPER_COHORT/'
basket_scores_path = Path(report_dir) / Path('basket_scores_4th_gen_zscored.tsv')
protein_scores_path = Path(report_dir) / Path('full_proteome_measures_z.tsv')
#protein_scores_path = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.06.22_AhS_PAPER_COHORT/full_proteome_measures_z.tsv'
meta_data_path = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Metadata_Papercohort_230801.xlsx'

instructions = pd.read_excel('/home/amir/Desktop/Annika_files/Figure 4F heatmap.xlsx')
instructions.set_index('Sample  name',inplace=True)

patients_list = instructions.index
basket_names =['EGFR','AXL','VEGFR','MET']
protein_names = ['PTEN', 'EGFR',	'MET',	'SMARCB1',	'CDKN2A-p16;CDKN2B','MEN1']
protein_scores_df = pd.read_csv(protein_scores_path,sep='\t')
protein_scores_df.columns = protein_scores_df.columns.str.replace('zscore_','')
protein_scores_df.set_index('Gene names',inplace=True)
protein_scores_df = protein_scores_df.loc[protein_names,patients_list]

basekt_scores_df = pd.read_csv(basket_scores_path,sep='\t')
basekt_scores_df.set_index('Sample',inplace=True)
basekt_scores_df = basekt_scores_df.loc[patients_list,basket_names].T


basekt_scores_df.index = basekt_scores_df.index + ' TUPAC_score'
protein_scores_df.index = protein_scores_df.index + ' Z_score'
final = pd.concat([basekt_scores_df,protein_scores_df],axis=0)



meta_data = pd.read_excel(meta_data_path)
meta =  meta_data[meta_data['Sample name'].isin(patients_list)][['Sample name', 'Paper_pseudo_identifier']].set_index('Sample name')
meta_dic = meta_data[meta_data['Sample name'].isin(patients_list)][['Sample name', 'Paper_pseudo_identifier']].set_index('Sample name').to_dict()['Paper_pseudo_identifier']


df = final.T.copy()
ps_names = meta.loc[df.index.tolist(),:]
#df.index = df.index.map(meta_dic)
min_df = df.min().min()
max_df = df.max().max()
print(max_df)
print(min_df)
df = df.fillna(-4.2)

plt.rcParams["figure.figsize"] = [20, 17]
plt.rcParams["figure.autolayout"] = True

cmap = colors.ListedColormap(['green','darkblue','lightblue','orange','red','darkred'])
bounds = [min_df,-2,0,1,1.5,2,max_df]
norm = colors.BoundaryNorm(bounds,cmap.N)

# plot heatmap
sns.heatmap(df,cmap=cmap,norm=norm)

ps_names.to_excel('/home/amir/Desktop/Annika_files/pseudo_names.xlsx')
plt.savefig('/home/amir/Desktop/Annika_files/figure4F_original_names_paper_cohort.svg')
plt.show()
