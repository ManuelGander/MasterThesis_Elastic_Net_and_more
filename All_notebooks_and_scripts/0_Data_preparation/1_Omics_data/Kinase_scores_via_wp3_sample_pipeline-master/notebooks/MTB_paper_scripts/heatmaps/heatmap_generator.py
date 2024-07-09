import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors


# example : python heatmap_generator.py  --t data.toml 
def make_the_scores_data_frame(toml_data):
    list_patients = toml_data['data']['list_patients']
    list_peptides = toml_data['data']['list_peptides']
    list_proteins = toml_data['data']['list_proteins']
    report_directory = toml_data['data']['report_directory']
    kinase_scores_df = pd.read_csv(Path(report_directory) /'kinase_results' / 'kinase_scores.tsv',sep='\t')
    phospho_scores = pd.read_csv(Path(report_directory) /'protein_results' / 'protein_scores.tsv',sep='\t')
    phospho_scores.set_index('Gene names',inplace=True)
    FP_scores = pd.read_csv(Path(report_directory) /'full_proteome_measures_z.tsv',sep='\t')
    phospho_peptides_scores = pd.read_csv(Path(report_directory) /'phospho_measures_z.tsv',sep='\t')
    FP_scores.set_index('Gene names',inplace=True)
    kinase_scores_df.set_index('PSP Kinases',inplace=True)
    phospho_peptides_scores.set_index('Modified sequence',inplace=True)
    FP_scores_columns = ['zscore_' + x for x in list_patients]
    FP_scores = FP_scores[FP_scores_columns]
    FP_scores = FP_scores[FP_scores.index.isin(list_proteins)]
    phospho_peptides_scores = phospho_peptides_scores[FP_scores_columns]
    phospho_peptides_scores = phospho_peptides_scores[phospho_peptides_scores.index.isin(list_peptides)]
    kinase_scores_df = kinase_scores_df[list_patients]
    kinase_scores_df = kinase_scores_df[kinase_scores_df.index.isin(list_proteins)]
    phospho_scores = phospho_scores[list_patients]
    phospho_scores = phospho_scores[phospho_scores.index.isin(list_proteins)]
    FP_scores.index = FP_scores.index + ' Z_score'
    phospho_peptides_scores.index = phospho_peptides_scores.index + ' Z_score'
    kinase_scores_df.index = kinase_scores_df.index + ' Kinase_Score'
    phospho_scores.index = phospho_scores.index + ' Phospho_Score'
    FP_scores.columns = FP_scores.columns.str.replace('zscore_','')
    phospho_peptides_scores.columns = phospho_peptides_scores.columns.str.replace('zscore_','')
    return pd.concat([FP_scores,phospho_peptides_scores,kinase_scores_df,phospho_scores])


def make_heatmap_plot(toml_data,df):
    report_directory = toml_data['data']['report_directory']
    min_df = df.min().min()
    max_df = df.max().max()
    grey_margin = min_df - 1
    abs_min =  grey_margin - 1
    df = df.fillna(abs_min)
    plt.rcParams["figure.figsize"] = [20, 7]
    plt.rcParams["figure.autolayout"] = True
    cmap = colors.ListedColormap(['grey','blue','red'])
    bounds = [abs_min,grey_margin,2,max_df]
    norm = colors.BoundaryNorm(bounds,cmap.N)
    # plot heatmap
    sns.heatmap(df,cmap=cmap,norm=norm)
    plt.savefig(Path(report_directory) / 'Heatmap_kinase_scores_phosopho_scores.svg')
    print(f'saved in {report_directory}/Heatmap_kinase_scores_phosopho_scores.svg')
    #plt.show()

if __name__ == '__main__':
    import sys
    import json
    import argparse
    import toml

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tomlfile", required=True,
                        help="Absolute path to configuration file.")

    args = parser.parse_args(sys.argv[1:])
    file_name = args.tomlfile
    toml_data=toml.load(file_name)
    print(toml_data)


    df = make_the_scores_data_frame(toml_data)
    make_heatmap_plot(toml_data,df)

  