import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', None)

#folder = '2022.06.30_test_miniml'
#folder = '2022.06.02_MT_summed_zscore_2nd_gen_baskets'
#folder = '2022.07.13_test'
#folder = '2022.07.13_test_MT'
#folder = '2022.07.13_test_MT_regular_median_centering'
#folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq'
#folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq_0.3missingness'

#folder = '2022.07.26_CJ_new_median_centering_and_lfq'
#folder = '2023.04.03_CJ_batch93_94'
folder = '2023.04.26_CJ_final_cohort2'
data_type = 'phospho'
#data_type = 'full_proteome'

df7 = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/{data_type}_measures_z.tsv', sep='\t')

median_vs_count_df = pd.DataFrame(pd.concat([df7.median(axis=0), df7.count()],axis=1))
median_vs_count_df = median_vs_count_df.rename(columns={0:'Median z-score', 1:'Number of IDs'})

print(median_vs_count_df)

median_vs_count_df.plot.scatter(x='Median z-score', y='Number of IDs')

plt.xlim([-1.25, 1.25])
plt.show()
