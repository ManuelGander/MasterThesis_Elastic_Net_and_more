import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', None)

#folder = '2022.06.30_test_miniml'
#folder = '2022.06.02_MT_summed_zscore_2nd_gen_baskets'
#folder = '2022.07.13_test'
#folder = '2022.07.13_test_MT'
#folder = '2022.07.13_test_MT_regular_median_centering'
#folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq'
folder_1 = '2022.08.02_MT_mixed_cohort'
folder_2 = '2022.08.02_MT_mixed_cohort_with_SIMSI'

df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_1}/preprocessed_fp.csv', index_col='Gene names')
df = df.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+')

df2 = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_2}/preprocessed_fp.csv', index_col='Gene names')
df2 = df2.filter(regex=r'^[A-Z,a-z]+\d{1,3}-\S+-\S+')

# to compare genes before and after SIMSI: 
# pd.concat([df.loc['PLA2G4C'], df2.loc['PLA2G4C']], axis=1)
# pd.concat([df.loc['PLA2G4C'].rename('without SIMSI'), df2.loc['PLA2G4C'].rename('with SIMSI')], axis=1).plot.scatter(x='without SIMSI', y='with SIMSI')

occ_before = df.count(axis=1).to_frame('occurrence_without_SIMSI')
median_before = df.median(axis=1).to_frame('median_intensity_without_SIMSI')
occ_after = df2.count(axis=1).to_frame('occurrence_with_SIMSI')
median_after = df2.median(axis=1).to_frame('median_intensity_with_SIMSI')

# print(df)

below_intensity_cutoff_before = (df < 3.5).sum(axis=1).to_frame('below_intensity_threshold_without_SIMSI')
below_intensity_cutoff_after = (df2 < 3.5).sum(axis=1).to_frame('below_intensity_threshold_with_SIMSI')

merged_df = pd.DataFrame().join([occ_before, occ_after, median_before, median_after, below_intensity_cutoff_before, below_intensity_cutoff_after], how='outer')

merged_df = merged_df.reset_index().rename(columns={'index': 'Gene names'})
merged_df.to_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_2}/protein_occurrences_and_median_abundances.tsv', sep='\t', index=False)

merged_df = merged_df.set_index('Gene names')
merged_df = merged_df.iloc[:, [0,1,2,3]]
merged_df1 = merged_df.iloc[:, [0,1]]
merged_df2 = merged_df.iloc[:, [2,3]]


fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(15,20))
fig.subplots_adjust(wspace=0.01)
p1 = sns.heatmap(merged_df1, cmap=plt.cm.get_cmap('Blues'), ax=ax1) #, cbar=False)
p1.figure.axes[-1].set_ylabel('Occurence in cohort', size=20)# cbaxes = fig.add_axes([0.01, 0.01, 0.03, 0.9])
p1.figure.axes[-1].tick_params(labelsize=14)
# fig.colorbar(ax1.collections[0], ax=ax1, reversed=True) #, use_gridspec=False, pad=0.5, cax=cbaxes)
p2 = sns.heatmap(merged_df2, cmap=plt.cm.get_cmap('rocket').reversed(), ax=ax2) #, cbar=False)
p2.figure.axes[-1].set_ylabel('Median intensity', size=20)
p2.figure.axes[-1].tick_params(labelsize=14)

# axp = ax1.imshow(np.random.randint(0, 300, (10, 10)))

# cb = plt.colorbar(axp, cax = cbaxes)![](../../../../kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.08.02_MT_mixed_cohort_with_SIMSI/protein_occurrences_and_median_abundances.jpg)
# fig.colorbar(ax2.collections[0], ax=ax2, location='right', use_gridspec=False, pad=0.5)
# ax2.yaxis.tick_right()
# ax1.set_yticklabels(labels=ax1.g.index, size=14)
# ax2.set_yticklabels(labels=merged_df2.index, size=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_ylabel('Gene names', size=16)
ax2.set_ylabel('')
ax1.set_xticklabels(labels=merged_df1.columns, rotation=45)
ax2.set_xticklabels(labels=merged_df2.columns, rotation=45)
plt.suptitle('Missingness vs intensities', fontsize=28, y=1)
plt.tight_layout()
plt.savefig(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder_2}/protein_occurrences_and_median_abundances.jpg', bbox_inches='tight', dpi=400)

