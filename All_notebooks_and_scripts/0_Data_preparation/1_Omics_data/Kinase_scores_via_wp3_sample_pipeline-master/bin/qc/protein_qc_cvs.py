import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

pd.set_option('expand_frame_repr', False)
pd.set_option('display.max_rows', None)


def main():
    #folder = '2022.06.30_test_miniml'
    #folder = '2022.06.02_MT_summed_zscore_2nd_gen_baskets'
    #folder = '2022.07.13_test'
    #folder = '2022.07.13_test_MT'
    #folder = '2022.07.13_test_MT_regular_median_centering'
    #folder = '2022.07.20_MT_mixed_cohort_new_median_centering_and_maxlfq'
    #folder = '2022.08.02_MT_mixed_cohort'
    #folder = '2022.08.02_MT_mixed_cohort_with_SIMSI'
    #folder = '2022.07.28_MT_add_single_peptide_id_info'
    folder = '2022.08.12_MT_inter_batch_CV'

    fp_df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_fp2.csv', index_col='Gene names')
    xlim = [5,10]
    plot_CVs(fp_df, xlim, min_observed=20)

    pp_df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_pp2.csv', index_col='Modified sequence')
    xlim = [4.5,8.5]
    plot_CVs(pp_df, xlim, min_observed=20)

    df2 = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/preprocessed_pp2_after_ms1_correction.csv', index_col='Modified sequence')
    df2 = df2.sort_index()
    df2 = convert_long_to_wide_format(create_metadata_columns(df2.reset_index()))
    xlim = [4.5,8.5]
    plot_CVs(df2, xlim, min_observed=20, input_is_log10=False)


def plot_CVs(df, xlim, min_observed=None, input_is_log10=True):
    #df = df.filter(regex='^Reporter intensity corrected (9|10|11)')
    #df = df.filter(regex='^Reporter intensity corrected')
    df = df.filter(regex='^Reporter intensity corrected (10|11)')
    if input_is_log10:
        df = np.power(10, df)
    
    #df = df.replace(0, np.nan)
    
    if min_observed is None or min_observed > df.shape[1]:
        df = df.dropna()
    else:
        df = df[df.filter(regex='^Reporter intensity corrected').count(axis=1) >= min_observed] # at least 10 channels are filled
    
    cv_before = df.std(axis=1) / df.mean(axis=1)
    
    #medians = df.median(axis=0)
    #median_correction = medians.mean() / medians
    #df = df.mul(median_correction, axis=1)
    
    cv = df.std(axis=1) / df.mean(axis=1)
    
    df = df.apply(np.log10)
    df['Median intensity'] = df.median(axis=1)
    df['CV'] = cv
    df = df[~df['Median intensity'].isna()]
    df = df[~df['CV'].isna()]
    
    print(df.head(n=100))
    print("Fraction with log10 intensity < 3.5: ", sum(df['Median intensity'] < 3.5) / len(df.index))
    
    #df.plot.scatter(x='Median intensity', y='CV')
    #plt.show()
    
    ax = density_scatter(np.array(df['Median intensity']), np.array(df['CV']), bins=100, s=1)
    plotDensityScatterTrendLine(np.array(df['Median intensity']), np.array(df['CV']), ax, numBins = 20)
    ax.set_xlabel('log10(intensity)')
    ax.set_ylabel('CV')
    
    ax.plot(xlim, [0.2, 0.2], '--')
    plt.tight_layout()
    plt.xlim(xlim)
    plt.ylim([0,1.5])
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



def plotDensityScatterTrendLine(xdat, ydat, ax, numBins = 50):
  x, y = zip(*[(x1,y1) for x1,y1 in sorted(zip(xdat, ydat))])
  binSize = len(x) / numBins
  means, stds, ints = list(), list(), list()
  for b in range(numBins):
    #mean, rmse = huber(y[b*binSize:(b+1)*binSize])
    means.append(np.median(y[int(b*binSize):int((b+1)*binSize)]))
    stds.append(np.std(y[int(b*binSize):int((b+1)*binSize)]))
    ints.append(x[int(np.ceil((b+0.5)*binSize))])
  ax.errorbar(ints, means, yerr = stds, linewidth=2.0, color = 'red')



if __name__ == "__main__":
    main()