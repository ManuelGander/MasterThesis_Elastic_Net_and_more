import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', None)

folder = '2022.07.26_MT_new_median_centering_and_lfq'

#data_type = 'phospho'
data_type = 'full_proteome'

df = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/pickedGeneGroups_with_quant.txt', sep='\t')


df = df.filter(regex='Unique peptides Reporter intensity corrected*')

num_proteins_df = pd.DataFrame(pd.concat([df.replace(0, np.nan).count(axis=0), df[df == 1].sum()],axis=1))
num_proteins_df = num_proteins_df.rename(columns={0:'Number of proteins', 1:'Number of proteins with 1 peptide'})

print(num_proteins_df)

num_proteins_df.plot.bar()

plt.tight_layout()
plt.show()
