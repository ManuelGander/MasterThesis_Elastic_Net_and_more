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

folder = '2022.07.26_CJ_new_median_centering_and_lfq'

df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2.csv')

batches = {re.sub(r'Reporter intensity corrected \d{1,2} ', '', c) for c in df.columns if c.startswith('Reporter intensity corrected')}

for batch in batches:
    tmt_cols = df.filter(regex=fr'^Reporter intensity corrected \d{{1,2}} {batch}$', axis='columns')
    df[f'Num channels detected {batch}'] = tmt_cols.count(axis='columns')
    imputation_cols = tmt_cols.columns.str.replace('Reporter intensity corrected', 'Identification metadata')
    df[imputation_cols] = np.where(tmt_cols.isna().multiply(df[f'Num channels detected {batch}'] > 0, axis='index'), 'n.d.', '')

df['Num batches detected'] = (df.filter(regex=r'^Num channels detected ') > 0).sum(axis='columns')
df['Num channels detected'] = df.filter(regex=r'^Num channels detected ').sum(axis='columns')
df['Num channels in batches detected'] = df['Num batches detected']*11
df['Fraction detected in batch but not in channel'] = 1.0 - df['Num channels detected'] / df['Num channels in batches detected']

df.sort_values(by='Fraction detected in batch but not in channel', ascending=False)[['Gene names', 'Num batches detected', 'Num channels detected', 'Num channels in batches detected', 'Fraction detected in batch but not in channel']].to_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/Fraction_detected_in_batch_but_not_in_channel.tsv', sep='\t', index=False)

df['Fraction detected in batch but not in channel'].replace(0.0, np.nan).plot.hist(bins=np.linspace(0,1,25))
plt.xlabel('Fraction detected in batch but not in channel')
plt.show()
