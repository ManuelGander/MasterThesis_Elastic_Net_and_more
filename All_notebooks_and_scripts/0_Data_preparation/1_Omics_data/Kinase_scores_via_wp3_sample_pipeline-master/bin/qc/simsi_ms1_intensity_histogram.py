import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

pval = 'p10'

folder = 'summaries_2022.08.02_MT_mixed_cohort_with_SIMSI'

df = pd.read_csv(
    f'/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI/simsi_output/FP/{folder}/p10/p10_evidence.txt',
    sep='\t',
    usecols=['Intensity', 'MS/MS count', 'Transferred spectra count'])

df = df[df['Intensity'] > 0]
# %% Transfers with quantified precursor

transferred_df = df[(df['Transferred spectra count'] > 0) & (df['Transferred spectra count'] == df['MS/MS count'])]
direct_df = df[(df['Transferred spectra count'] == 0)]

print(len(transferred_df.index))
print(len(direct_df.index))

bins = 100
alpha = 0.6
density = True


plt.hist(np.log10(direct_df['Intensity']), bins=bins, density=density, label='direct', alpha=alpha)
plt.hist(np.log10(transferred_df['Intensity']), bins=bins, density=density, label='transferred', alpha=alpha)
if density:    
    plt.ylim(0, 1)
plt.xlim(4, 11)
plt.xlabel('Log10 Intensity')
plt.ylabel('Count')
plt.legend()
plt.tight_layout()
plt.show()
