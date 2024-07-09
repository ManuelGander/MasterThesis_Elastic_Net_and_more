import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', None)

folder = '2022.07.13_test_MT'

#df4 = pd.read_csv('/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_pp2_before_ms1_correction.csv')
#df4 = pd.read_csv('/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_pp2.csv')
#df4 = pd.read_csv('/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2_before_ms1_correction.csv')
df4 = pd.read_csv('/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp.csv')

df4[df4['Batch'] == 'Chordoma_Batch2'][df4.columns[df4.columns.str.startswith('Reporter intensity corrected')]].apply(np.log10).plot.hist(bins=100, alpha=0.3,histtype='step')

# z-scoring
# x.sub( x.median(axis=1), axis=0).divide(x.std(axis=1), axis=0)

plt.show()
