import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


pd.set_option('display.max_rows', None)

#folder = '2022.07.13_test_MT'
#folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq'
folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq_0.3missingness'

df4 = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2_after_ms1_correction.csv')

ax = df4[df4['Batch'] == 'Sarcoma_Batch1']['MS1'].apply(np.log10).plot.hist(ax=ax, bins=100, alpha=0.3, histtype='step')
df4[df4['Batch'] == 'Chordoma_Batch1']['MS1'].apply(np.log10).plot.hist(bins=100, alpha=0.3, histtype='step')
df4[df4['Batch'] == 'Chordoma_Batch2']['MS1'].apply(np.log10).plot.hist(ax=ax, bins=100, alpha=0.3, histtype='step')

#>>> df4[df4['Batch'] == 'Sarcoma_Batch1']['MS1'].count()
#136990
#>>> df4[df4['Batch'] == 'Chordoma_Batch1']['MS1'].count()
#99048
#>>> df4[df4['Batch'] == 'Chordoma_Batch2']['MS1'].count()
#115269

# z-scoring
# z = x.sub( x.median(axis=1), axis=0).divide(x.std(axis=1), axis=0)

#Reporter intensity corrected 1 Sarcoma_Batch1      0.230332
#Reporter intensity corrected 2 Sarcoma_Batch1      0.304629
#Reporter intensity corrected 3 Sarcoma_Batch1      0.279266
#Reporter intensity corrected 4 Sarcoma_Batch1      0.233614
#Reporter intensity corrected 5 Sarcoma_Batch1      0.345547
#Reporter intensity corrected 6 Sarcoma_Batch1      0.244645
#Reporter intensity corrected 7 Sarcoma_Batch1      0.421922
#Reporter intensity corrected 8 Sarcoma_Batch1      0.264761
#Reporter intensity corrected 9 Sarcoma_Batch1      0.547925
#Reporter intensity corrected 10 Sarcoma_Batch1     0.595746
#Reporter intensity corrected 11 Sarcoma_Batch1     0.536733

#Reporter intensity corrected 1 Chordoma_Batch1    -0.390634
#Reporter intensity corrected 2 Chordoma_Batch1    -0.395586
#Reporter intensity corrected 3 Chordoma_Batch1    -0.381363
#Reporter intensity corrected 4 Chordoma_Batch1    -0.443216
#Reporter intensity corrected 5 Chordoma_Batch1    -0.394203
#Reporter intensity corrected 6 Chordoma_Batch1    -0.261735
#Reporter intensity corrected 7 Chordoma_Batch1    -0.374128
#Reporter intensity corrected 8 Chordoma_Batch1    -0.332365
#Reporter intensity corrected 9 Chordoma_Batch1    -0.077152
#Reporter intensity corrected 10 Chordoma_Batch1   -0.030257
#Reporter intensity corrected 11 Chordoma_Batch1   -0.157961

#Reporter intensity corrected 1 Chordoma_Batch2    -0.226986
#Reporter intensity corrected 2 Chordoma_Batch2    -1.090985
#Reporter intensity corrected 3 Chordoma_Batch2    -0.192034
#Reporter intensity corrected 4 Chordoma_Batch2    -0.038891
#Reporter intensity corrected 5 Chordoma_Batch2    -0.099182
#Reporter intensity corrected 6 Chordoma_Batch2    -0.081537
#Reporter intensity corrected 7 Chordoma_Batch2     0.004823
#Reporter intensity corrected 8 Chordoma_Batch2    -0.176301
#Reporter intensity corrected 9 Chordoma_Batch2     0.138920
#Reporter intensity corrected 10 Chordoma_Batch2    0.230607
#Reporter intensity corrected 11 Chordoma_Batch2    0.082834

# Sarcoma_Batch1
#Reporter intensity corrected 1     130349
#Reporter intensity corrected 2     134307
#Reporter intensity corrected 3     134458
#Reporter intensity corrected 4     134533
#Reporter intensity corrected 5     136164
#Reporter intensity corrected 6     134288
#Reporter intensity corrected 7     134449
#Reporter intensity corrected 8     131893
#Reporter intensity corrected 9     134304
#Reporter intensity corrected 10    128339
#Reporter intensity corrected 11    126161

# Chordoma_Batch1
#Reporter intensity corrected 1      94419
#Reporter intensity corrected 2      96503
#Reporter intensity corrected 3      96180
#Reporter intensity corrected 4      94279
#Reporter intensity corrected 5      97062
#Reporter intensity corrected 6      97325
#Reporter intensity corrected 7      97797
#Reporter intensity corrected 8      97701
#Reporter intensity corrected 9      92339
#Reporter intensity corrected 10     90353
#Reporter intensity corrected 11     87035

# Chordoma_Batch2
#Reporter intensity corrected 1     105382
#Reporter intensity corrected 2      76723
#Reporter intensity corrected 3     111832
#Reporter intensity corrected 4     111130
#Reporter intensity corrected 5     110648
#Reporter intensity corrected 6     110756
#Reporter intensity corrected 7     114538
#Reporter intensity corrected 8     112903
#Reporter intensity corrected 9     106960
#Reporter intensity corrected 10    107640
#Reporter intensity corrected 11    102558


def drop_duplicate_indices(df):
    return df[~df.index.duplicated(keep='first')]

key = ['Modified sequence', 'Charge']
#sarcoma1_df = df4[df4['Batch'] == 'Sarcoma_Batch20'].set_index(key)
#sarcoma1_df = drop_duplicate_indices(sarcoma1_df)

#chordoma1_df = df4[df4['Batch'] == 'Chordoma_Batch1'].set_index(key)
#chordoma1_df = drop_duplicate_indices(chordoma1_df)

#chordoma2_df = df4[df4['Batch'] == 'Chordoma_Batch2'].set_index(key)
#chordoma2_df = drop_duplicate_indices(chordoma2_df)

#merged_df = sarcoma1_df[['MS1']].join(chordoma1_df[['MS1']], lsuffix="_sarcoma1", rsuffix="_chordoma1", how='outer').join(chordoma2_df[['MS1']], how='outer').rename(columns={'MS1': 'MS1_chordoma2'})
#merged_df = merged_df.replace(0.0, np.nan).apply(np.log10)

#intensity_cols = ['Reporter intensity corrected 2']
intensity_cols = [x for x in df4.columns if x.startswith('Reporter intensity corrected')]
#intensity_cols = 'MS1'

ms1s = df4.set_index(key)[intensity_cols + ['Batch']]
all_batches = {}
for batch_name, df in ms1s.groupby('Batch'):
    df = drop_duplicate_indices(df)
    df = df.drop(columns='Batch')
    df = df.rename(columns=lambda x: f'{x} {batch_name}')
    all_batches[batch_name] = df

merged_df = pd.DataFrame().join(all_batches.values(), how="outer")
merged_df = merged_df.replace(0.0, np.nan).apply(np.log10)

merged_df_filtered = merged_df[merged_df.count(axis=1) > 0.7*len(merged_df.columns)]
#diff = merged_df['MS1_chordoma1'] - merged_df['MS1_chordoma2']
#diff.median()
#plt.hist(diff, bins=np.linspace(-2,2,100))


merged_df.plot.hist(histtype='step', bins=30)

merged_df.dropna().plot.hist(histtype='step', bins=30)

ax = merged_df.dropna()['MS1_chordoma1'].plot.hist(histtype='step', bins=30, label='all_sarcoma1')
merged_df[merged_df['MS1_chordoma1'].isna() | merged_df['MS1_chordoma2'].isna()]['MS1_chordoma1'].plot.hist(histtype='step', bins=30, ax=ax, label='exclusive_sarcoma1')
plt.legend()
plt.show()

