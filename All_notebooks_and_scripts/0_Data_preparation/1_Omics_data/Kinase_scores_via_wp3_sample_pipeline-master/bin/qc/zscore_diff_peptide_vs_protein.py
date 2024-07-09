import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

pd.set_option('display.max_rows', None)


from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn


def main(argv):
    #folder = '2022.07.13_test_MT'
    folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq_0.3missingness'

    peptide_abundances = get_peptide_abundances(folder)
    peptide_zscores = get_zscores(peptide_abundances)
    #plot_median_vs_count(peptide_zscores)

    protein_abundances = get_protein_abundances(folder)
    protein_zscores = get_zscores(protein_abundances)
    #plot_median_vs_count(protein_zscores)

    protein_counts = get_protein_abundances(folder, col_start='Reporter intensity count').apply(lambda x : np.power(10, x))

    peptide_protein_merged_df = protein_zscores.join(peptide_zscores, lsuffix="_protein_z", rsuffix="_peptide_z", how='inner').join(protein_counts)
    #peptide_protein_merged_df = peptide_protein_merged_df.dropna(axis=0)
        
    tmt_channel = 2
    #batch = 'Chordoma_Batch1'
    #batch = 'Chordoma_Batch2'
    #batch = 'Sarcoma_Batch1'
    batch = 'Sarcoma_Batch30'

    plot_peptide_protein_diff_histogram(peptide_protein_merged_df, tmt_channel, batch)


    tmt_channel_1 = 3
    tmt_channel_2 = 3
    batch_1 = 'Sarcoma_Batch20'
    batch_2 = 'Sarcoma_Batch3'
    #batch = 'Sarcoma_Batch1'

    #for batch_2 in ['Sarcoma_Batch2', 'Sarcoma_Batch3', 'Sarcoma_Batch4', 'Sarcoma_Batch23', 'Sarcoma_Batch24', 'Sarcoma_Batch25', 'Sarcoma_Batch30', 'Sarcoma_Batch31', 'Sarcoma_Batch32']:
    plot_peptide_abundance_correlation(peptide_abundances, tmt_channel_1, batch_1, tmt_channel_2, batch_2)
    plot_peptide_protein_diff_vs_count_scatter(peptide_protein_merged_df, tmt_channel_1, batch_1, tmt_channel_2, batch_2)


def plot_median_vs_count(df7):
    median_vs_count_df = pd.DataFrame(pd.concat([df7.median(axis=0), df7.count()],axis=1))
    median_vs_count_df = median_vs_count_df.rename(columns={0:'Median z-score', 1:'Number of IDs'})
    
    print(median_vs_count_df)
    
    median_vs_count_df.plot.scatter(x='Median z-score', y='Number of IDs')
    
    plt.xlim([-1.25, 1.25])
    plt.show()


def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)
    
    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0
    
    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    
    ax.scatter( x, y, c=z, **kwargs )
    
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Density')
    
    return ax
    
def plot_peptide_protein_diff_histogram(peptide_protein_merged_df, tmt_channel, batch):
    print(batch, tmt_channel)
    diff = peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel} {batch}_protein_z'] - peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel} {batch}_peptide_z']
    
    print(diff.mean())
    plt.hist(diff, bins=np.linspace(-3, 3, 50))
    plt.show()


def plot_peptide_abundance_correlation(peptide_abundances, tmt_channel_1, batch_1, tmt_channel_2, batch_2):
    merged_peptide_df = pd.DataFrame({'peptide_1': peptide_abundances[f'Reporter intensity corrected {tmt_channel_1} {batch_1}'], 'peptide_2': peptide_abundances[f'Reporter intensity corrected {tmt_channel_2} {batch_2}']})
    
    bins = np.linspace(2, 10, 50)
    ax = merged_peptide_df[merged_peptide_df['peptide_1'].isna()]['peptide_2'].plot.hist(histtype='step', label=f'unique {tmt_channel_2} {batch_2}', bins=bins)
    merged_peptide_df[merged_peptide_df['peptide_2'].isna()]['peptide_1'].plot.hist(ax = ax, histtype='step', label=f'unique {tmt_channel_1} {batch_1}', bins=bins)
    merged_peptide_df.dropna()['peptide_1'].plot.hist(ax = ax, histtype='step', label=f'common {tmt_channel_1} {batch_1}', bins=bins)
    merged_peptide_df.dropna()['peptide_2'].plot.hist(ax = ax, histtype='step', label=f'common {tmt_channel_2} {batch_2}', bins=bins)
    ax.legend()
    
    merged_peptide_df = merged_peptide_df.dropna()
    ax = density_scatter(merged_peptide_df['peptide_1'], merged_peptide_df['peptide_2'], s=1)
    ax.plot([2, 10], [2,10], 'r-')
    ax.set_xlim([2,10])
    ax.set_ylim([2,10])
    ax.set_xlabel(f'Reporter intensity corrected {tmt_channel_1} {batch_1}')
    ax.set_ylabel(f'Reporter intensity corrected {tmt_channel_2} {batch_2}')
    plt.tight_layout()
    plt.show()
    #plt.savefig(f'/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/2022.07.14_MT_low_ids_low_basket_scores_investigation/peptide_correlation_batch_{batch_1}_vs_{batch_2}.png')


def plot_peptide_protein_diff_vs_count_scatter(peptide_protein_merged_df, tmt_channel_1, batch_1, tmt_channel_2, batch_2):
    print(batch_1, tmt_channel_1)
    diff = peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel_1} {batch_1}_protein_z'] - peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel_2} {batch_2}_protein_z']
    
    diff_peptide = peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel_1} {batch_1}_peptide_z'] - peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel_2} {batch_2}_peptide_z']
    
    diff_count = peptide_protein_merged_df[f'Reporter intensity count {tmt_channel_1} {batch_1}'] - peptide_protein_merged_df[f'Reporter intensity count {tmt_channel_2} {batch_2}']
    
    diff_z_similar_counts = diff[np.abs(diff_count) < 3]
    
    print(diff.dropna().count())
    print(diff.mean())
    print(diff_z_similar_counts.mean())
    print(diff_count.mean())
    
    merged_diff_df = pd.DataFrame({'diff': diff, 'diff_peptide': diff_peptide, 'diff_count': diff_count})
    merged_diff_df = merged_diff_df.drop_duplicates().dropna()
    
    density_scatter(merged_diff_df['diff'], merged_diff_df['diff_count'], s=1)
    
    plt.xlabel("z-score diff")
    plt.ylabel("peptide count diff")
    
    max_zscore_plot = 5
    max_count_diff_plot = 50
    plt.xlim([-1*max_zscore_plot, max_zscore_plot])
    plt.ylim([-1*max_count_diff_plot, max_count_diff_plot])
    
    merged_diff_agg_df = merged_diff_df.groupby(level=0).agg('mean')
    ax = density_scatter(merged_diff_agg_df['diff'], merged_diff_agg_df['diff_peptide'], s=2)
    ax.plot([-5, 5], [-5,5], 'r-')
    plt.xlim([-1*max_zscore_plot, max_zscore_plot])
    plt.ylim([-1*max_zscore_plot, max_zscore_plot])
    
    merged_peptide_df = pd.DataFrame({'peptide_1': peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel_1} {batch_1}_peptide_z'], 'peptide_2': peptide_protein_merged_df[f'Reporter intensity corrected {tmt_channel_2} {batch_2}_peptide_z']})
merged_peptide_df = merged_peptide_df.dropna()
density_scatter(merged_peptide_df['peptide_1'], merged_peptide_df['peptide_2'], s=1)
    
    plt.figure()
    plt.hist(diff, bins=30, histtype='step', label='All proteins')
    plt.hist(diff_z_similar_counts, bins=30, histtype='step', label='Similar peptide count')
    plt.xlim([-max_zscore_plot, max_zscore_plot])
    plt.legend()
    
    plt.figure()
    plt.hist(diff_count, bins=np.linspace(-1*max_count_diff_plot, max_count_diff_plot, 30))
    plt.xlim([-1*max_count_diff_plot, max_count_diff_plot])
    
    plt.show()


def get_peptide_abundances(folder, index_col='Gene names'):
    #df4 = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2_before_ms1_correction.csv')
    df4 = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2_after_ms1_correction.csv')
    if index_col == 'Gene names':
        df4[index_col] = df4['Gene names'].fillna("") + "_" + df4['Modified sequence']
    
    intensity_cols = [x for x in df4.columns if x.startswith('Reporter intensity corrected')]
    
    ms1s = df4.set_index(index_col)[intensity_cols + ['Batch']]
    all_batches = {}
    for batch_name, df in ms1s.groupby('Batch'):
        df = drop_duplicate_indices(df)
        df = df.drop(columns='Batch')
        df = df.rename(columns=lambda x: f'{x} {batch_name}')
        all_batches[batch_name] = df
    
    merged_df = pd.DataFrame().join(all_batches.values(), how="outer")
    merged_df = merged_df.replace(0.0, np.nan).apply(np.log10)
    
    if index_col == 'Gene names':
        merged_df.index = map(remove_mod_sequence, merged_df.index)
    else:
        merged_df.index = map(remove_modifications, merged_df.index)
    
    return merged_df


def get_protein_abundances(folder, index_col='Gene names', col_start='Reporter intensity corrected'):
    df = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2.csv')
    #df = df.set_index('Best peptide')
    df = df.set_index(index_col)
    
    df = drop_duplicate_indices(df)
    
    df = df.rename(columns=lambda x: x.replace('Unique peptides Reporter intensity corrected', 'Reporter intensity count'))
    
    tmt_cols_merged = [x for x in df.columns if x.startswith(col_start) and not '9' in x and not '10' in x and not '11' in x]
    return df[tmt_cols_merged].apply(np.log10)
    

def remove_modifications(mod_sequence, remove_phospho_only=False):
    raw_sequence = mod_sequence.replace('(Phospho (STY))', '')
    if remove_phospho_only:
        return raw_sequence
    raw_sequence = raw_sequence.replace('(Acetyl (Protein N-term))', '')
    raw_sequence = raw_sequence.replace('(Oxidation (M))', '')
    raw_sequence = raw_sequence.replace('_', '')
    return raw_sequence


def drop_duplicate_indices(df):
    return df[~df.index.duplicated(keep='first')]


def remove_mod_sequence(key):
    return key.split("_")[0]


def get_zscores(df):
    return df.sub(df.median(axis=1), axis=0).divide(df.std(axis=1), axis=0)


if __name__ == "__main__":
    main(sys.argv[1:])
