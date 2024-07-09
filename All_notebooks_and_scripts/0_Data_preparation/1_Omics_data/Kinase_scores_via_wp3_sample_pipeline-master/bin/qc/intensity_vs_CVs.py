import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


def main(argv):
    #batch = 'Batch1_FP_INFORM_MASTER' # Lumos 3
    batch = 'Batch39_FP_MASTER' # Lumos 3
    #batch = 'Batch46_FP_MASTER' # Lumos 1
    
    #batch = 'Batch1_PP_INFORM_MASTER' # Lumos 2
    #batch = 'Batch39_PP_MASTER' # Lumos 2
    #batch = 'Batch46_PP_MASTER' # Lumos 2
    
    msms_df = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Searches/Sarcoma/{batch}/combined/txt/msmsScans.txt', sep='\t')
    xlim = [2.5,6.5]
    plot_CVs(msms_df, xlim)
        
    proteins_df = pd.read_csv(f'/home/matthewt/TOPAS/WP31/Searches/Sarcoma/{batch}/combined/txt/proteinGroups.txt', sep='\t')
    xlim = [2.5,7.5]
    plot_CVs(proteins_df, xlim)


def plot_CVs(df, xlim):
    #df = df.filter(regex='^Reporter intensity corrected (10|11)')
    df = df.filter(regex='^Reporter intensity corrected (9|10|11)')
    #df = msms_df.filter(regex='^Reporter intensity corrected')
    df = df.replace(0, np.nan)
    df = df.dropna()
    #df = df[df.filter(regex='^Reporter intensity corrected').count(axis=1) >= 11] # at least 5 out of 11 channels are filled
    
    cv_before = df.std(axis=1) / df.mean(axis=1)
    
    medians = df.median(axis=0)
    median_correction = medians.mean() / medians
    df = df.mul(median_correction, axis=1)
    
    cv = df.std(axis=1) / df.mean(axis=1)
    
    df = df.apply(np.log10)
    df['Median intensity'] = df.median(axis=1)
    df['CV'] = cv
    df = df[~df['Median intensity'].isna()]
    df = df[~df['CV'].isna()]
    
    print("Fraction with log10 intensity < 3.5: ", sum(df['Median intensity'] < 3.5) / len(df.index))
    
    #df.plot.scatter(x='Median intensity', y='CV')
    #plt.show()
    
    ax = density_scatter(np.array(df['Median intensity']), np.array(df['CV']), bins=50, s=1)
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
    means.append(np.mean(y[int(b*binSize):int((b+1)*binSize)]))
    stds.append(np.std(y[int(b*binSize):int((b+1)*binSize)]))
    ints.append(x[int(np.ceil((b+0.5)*binSize))])
  ax.errorbar(ints, means, yerr = stds, linewidth=2.0, color = 'red')


if __name__ == "__main__":
    main(sys.argv[1:])
