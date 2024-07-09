import sys
import logging

import pandas as pd

import bin.config as config
from bin.preprocess_tools import rename_columns_with_sample_ids
from bin.sample_annotation import get_channel_to_sample_id_dict

logger = logging.getLogger(__name__)

pd.set_option('display.max_rows', None)


def save_as_patient_isoform_table(results_folder, sample_annotation, gene_name):
    df = pd.read_csv(f"{results_folder}/proteinGroups_{gene_name}.txt", sep='\t')
    
    sample_annotation_df = pd.read_csv(sample_annotation)
    channel_to_sample_id_dict = get_channel_to_sample_id_dict(sample_annotation_df, remove_qc_failed=True)
    
    index_cols = ["Protein IDs", "Peptide counts (unique)"]
    df = rename_columns_with_sample_ids(df, channel_to_sample_id_dict, index_cols=index_cols)
    
    output_file = f"{results_folder}/proteinGroups_{gene_name}_patient_ids.txt"
    df.T.to_csv(output_file, sep='\t', header=False)
    
    print(output_file)
    

# python -m bin.qc.picked_group_fdr_single_gene_postprocessing ../../config_patients.json CDKN2A
if __name__ == "__main__":
    config_file = sys.argv[1]
    gene_name = sys.argv[2]
    configs = config.load(config_file)

    save_as_patient_isoform_table(configs["results_folder"], configs['sample_annotation'], gene_name)
