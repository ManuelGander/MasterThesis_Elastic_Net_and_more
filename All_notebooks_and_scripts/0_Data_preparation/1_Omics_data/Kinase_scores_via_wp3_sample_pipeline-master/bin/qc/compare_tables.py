import pandas as pd

results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/'

#folder1 = '2023.05.30_MT_minimal_test_remove_proteins_without_peptides'
folder1 = '2023.05.30_MT_minimal_test_PSP_kinases'
folder2 = '2023.05.30_MT_minimal_test_PSP_kinases_Uniprot_IDs'

#file_name = 'basket_scores_4th_gen.tsv'
file_name = 'subbasket_scores_RET.tsv'

df1 = pd.read_csv(f'{results_folder}/{folder1}/{file_name}', sep='\t')
df2 = pd.read_csv(f'{results_folder}/{folder2}/{file_name}', sep='\t')

print(df1.round(5).compare(df2.round(5)))
