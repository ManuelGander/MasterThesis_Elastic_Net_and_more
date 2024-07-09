import re

import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

folder = '2022.08.02_MT_mixed_cohort'
#folder = '2022.08.02_MT_mixed_cohort_with_SIMSI'

df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/annot_pp.csv', index_col='Modified sequence')

pY_df = df[~df['Site positions identified (MQ)'].isna() & (df['Site positions identified (MQ)'].str.contains('_Y'))]

pY_df['Occurrence'] = pY_df.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+').count(axis=1)

num_pYs_per_patient_df = pY_df.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+').count(axis=0)

pY_df['pY sites'] = pY_df['Site positions identified (MQ)'].apply(lambda x: ";".join([s for s in x.split(";") if '_Y' in s]))
#pY_df = pY_df.drop(columns=[c for c in pY_df.columns if c.startswith('Identification metadata') or re.match(r'^[A-Z,a-z]+\d{1,3}-\S+-\S+', c) or c == 'Unnamed: 0'])

pY_df = pY_df.reset_index()

cols = ['pY sites', 'Occurrence', 'Gene names', 'PSP Kinases', 'PSP_LT_LIT', 'PSP_MS_LIT', 'PSP_MS_CST', 'PSP_ON_FUNCTION', 'PSP_ON_PROCESS', 'PSP_ON_PROT_INTERACT', 'PSP_ON_OTHER_INTERACT', 'PSP_NOTES', 'PSP_URL', 'basket', 'rtk', 'Site positions identified (MQ)', 'Modified sequence', 'Proteins']

pY_df[cols].to_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/pY_sites.tsv', sep='\t', index=False)

num_pYs_per_patient_df.to_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/pY_sites_per_patient.tsv', sep='\t', index=True)


