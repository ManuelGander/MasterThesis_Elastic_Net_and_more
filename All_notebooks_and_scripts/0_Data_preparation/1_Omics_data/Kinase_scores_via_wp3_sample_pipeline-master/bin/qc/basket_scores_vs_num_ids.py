import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

pd.set_option('display.max_rows', None)

#folder = '2022.06.30_test_miniml'
#folder = '2022.06.02_MT_summed_zscore_2nd_gen_baskets'
#folder = '2022.07.13_test'
#folder = '2022.07.13_test_MT'
#folder = '2022.07.13_test_MT_regular_median_centering'
#folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq'
#folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq_0.3missingness'
#folder = '2022.08.02_MT_mixed_cohort'
folder = '2022.08.02_MT_mixed_cohort_with_SIMSI'

#data_type = 'FP - Scores'
#data_type = 'FP - RTK Scores'
#data_type = 'PP - Scores'
#data_type = 'PP - RTK Scores'
data_type = 'PP - TOPAS Drug Scores'

df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/basket_scores.tsv', sep='\t', index_col='Sample')
df2 = df.filter(regex=f'^{data_type}')

#for col in df.columns:
#    df.plot.scatter(x=col, y=f'{data_type}.num_identified')
#    plt.show()

cols = df2.columns.values.tolist()
n = len(cols)
heatmap = np.zeros((n, n))
for i, col in enumerate(cols):
    for j, col2 in enumerate(cols[i+1:]):
        r2, p_value = stats.pearsonr(df2[col], df2[col2])
        print(col, col2, r2, p_value)
        heatmap[i, i+j+1] = r2
        #if r2 > 0.6:
        #    df2.plot.scatter(x=col, y=col2)
        #    plt.show()

cols_short = [c.replace(f'{data_type}.', '') for c in cols]

plt.imshow(heatmap, cmap='bwr', interpolation='nearest', vmin=-1, vmax=1)
cbar = plt.colorbar()
cbar.set_label('Pearson correlation', rotation=270)
plt.xticks(np.arange(0, n, 1), cols_short, rotation=90)
plt.yticks(np.arange(0, n, 1), cols_short)
plt.tight_layout()
plt.show()
