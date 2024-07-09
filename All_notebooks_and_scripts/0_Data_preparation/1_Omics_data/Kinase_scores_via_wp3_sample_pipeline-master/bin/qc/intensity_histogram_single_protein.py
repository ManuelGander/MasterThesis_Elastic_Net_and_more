import pandas as pd
import os
import matplotlib.pyplot as plt

results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.12.06_AhS_all_mixed_cohort/'
gene_name = 'PTEN'

df = pd.read_csv(os.path.join(results_folder, 'preprocessed_fp.csv'), index_col='Gene names')

df.filter(regex=r'^[A-Z]\d{3}-\S+-\S+').loc[gene_name].plot.hist(bins=30)

plt.xlabel('log10(intensity)')
plt.show()
