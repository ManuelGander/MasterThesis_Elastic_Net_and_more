import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.11.11_maxpep_1/'

df = pd.read_csv(os.path.join(results_folder, 'preprocessed_fp.csv'))
df = df.set_index('Gene names')
df = df.drop(df.loc[:, df.columns.str.startswith('Identification')].columns, axis=1)

df_z = pd.read_csv(os.path.join(results_folder, 'full_proteome_measures_z.tsv'), sep='\t')
df_z = df_z.set_index('Gene names')
df_z = df_z[[col for col in df_z if col.startswith('zscore')]]
df_z.columns = df_z.columns.str.removeprefix('zscore_')

pois = ['CDK4', 'PTEN']

# subset to protein of interest
df = df.loc[pois, :]
df_z = df_z.loc[pois, :]
print(df.head())
print(df_z.head())

# plot
for poi in pois:
    x, y, y2 = df.columns.tolist(), df.loc[poi, :].values.tolist(), df_z.loc[poi, :].values.tolist()
    fig, ax = plt.subplots()
    ax.scatter(x,
               y,
               color="red",
               marker="o",
               s=10)

    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    # set y-axis label
    ax.set_ylabel("Intensities", color="red", fontsize=14)

    # twin object for two different y-axis on the sample plot
    ax2 = ax.twinx()
    # make a plot with different y-axis using second axis object
    ax2.scatter(x,
                y2,
                color="blue",
                marker="o",
                s=10)
    ax2.set_ylabel("Z-scores", color="blue", fontsize=14)
    # save the plot as a file
    fig.savefig(os.path.join(results_folder, 'Results_investigation_plots', f'{poi}_intensity_and_zscore_cohort.jpg'),
                format='jpeg',
                dpi=300,
                bbox_inches='tight')
    plt.close()
