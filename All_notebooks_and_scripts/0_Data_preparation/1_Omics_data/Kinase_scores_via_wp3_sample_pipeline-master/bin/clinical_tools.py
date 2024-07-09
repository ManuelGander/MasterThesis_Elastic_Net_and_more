import json
import logging
from itertools import compress
from pathlib import Path
from typing import List, Dict, Tuple, Union

import pandas as pd
import numpy as np

import psite_annotation as pa
from . import utils

logger = logging.getLogger(__name__)


def phospho_annot(df: pd.DataFrame,
                  pspFastaFile: Union[str, Path],
                  pspKinaseSubstrateFile: Union[str, Path],
                  pspAnnotationFile: Union[str, Path],
                  pspRegulatoryFile: Union[str, Path]) -> pd.DataFrame:
    """
    Phospho-site annotation of experimental data using in-house developed tool (MT) based mainly on Phosphosite Plus but also in vitro
    experiments from the lab of Ishihama.

    :param df: dataframe with measured peptide intensities
    :param pspFastaFile: file used for adding peptide and psite positions
    :param pspKinaseSubstrateFile: file used for adding kinase substrate annotation
    :param pspAnnotationFile: file used for adding annotations from Phosphosite Plus
    :param pspRegulatoryFile: file used for adding regulatory information
    """
    logger.info('Phosphosite annotation')

    logger.info(f'Phospho data before: {df.shape}')
    df.reset_index(level='Modified sequence', inplace=True)
    df = pa.addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True, returnAllPotentialSites=False)
    df = pa.addPSPKinaseSubstrateAnnotations(df, pspKinaseSubstrateFile, gene_name=True)
    df = pa.addPSPAnnotations(df, pspAnnotationFile)
    df = pa.addPSPRegulatoryAnnotations(df, pspRegulatoryFile)

    df['PSP_LT_LIT'] = df['PSP_LT_LIT'].apply(lambda x: max(x.split(';')))
    df['PSP_MS_LIT'] = df['PSP_MS_LIT'].apply(lambda x: max(x.split(';')))
    df['PSP_MS_CST'] = df['PSP_MS_CST'].apply(lambda x: max(x.split(';')))
    df.rename(columns={'Site positions': 'Site positions identified (MQ)'}, inplace=True)
    df = pa.addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True, returnAllPotentialSites=True)
    logger.info(f'Phospho data after: {df.shape}')

    return df


def add_psp_urls(pp: pd.DataFrame) -> pd.DataFrame:
    """
    Function to add column to dataframe with URLs to proteins/isoforms that the
    peptides in the data belongs to:  https://www.phosphosite.org/. It uses already
    annotated information from PSP to check if any annotation exists. If it does, the URL
    is created from template and concatenated with the uniprot ID.

    :param pp: df to annotate with URL to PhosphoSitePlus
    :return: df with added annotation of URL to PhosphoSitePlus
    """
    pp[['PSP_URL', 'PSP_URL_extra']] = pp[['Start positions', 'Proteins']].apply(add_url_columns, axis=1, result_type="expand")
    return pp


def add_url_columns(row) -> Tuple[str, List]:
    start_positions, proteins = row
    # create boolean list for (p-site, protein) pairs found in PSP or not
    # check for any modified peptide with start position different from -1
    found_psites = [int(value) > 0 for value in start_positions.split(';') if value != '']
    # row[0] is integer index of row and row[1] is column value
    proteins = list(compress(proteins.split(';'),
                             found_psites))

    URLs = list()
    main_url = ""
    main_found = False

    # There can be found more than one main protein URL but then the first is used
    # and the rest goes with the isoform URLs
    url_start = "https://www.phosphosite.org/uniprotAccAction?id="
    for index, protein in enumerate(proteins):
        # don't allow isoforms (recognizable by "-" in their protein IDs) as main protein
        if '-' not in protein and not main_found:
            main_url = "=HYPERLINK(\"" + str(url_start) + str(protein) + "\")"
            main_found = True
        else:
            URLs.append(str(url_start) + str(protein))

    return main_url, URLs


def prot_basket_annotation(df: pd.DataFrame,
                           prot_baskets: Union[str, Path],
                           data_type: str,
                           basket_type: str) -> Tuple[pd.DataFrame, Dict]:
    """
    Adds columns with basket annotations and weights to a dataframe

    :param df: dataframe with a 'Gene names' column to be annotated
    :param prot_baskets: list of path(s) to file(s) with annotations
    :param data_type: either 'pp' for phospho or 'fp' for full proteome
    :param basket_type: either 'basket', 'other', 'rtk' corresponding to the sheets in the annotation file
    """
    logger.info(f'Proteomics baskets annotation {data_type} {basket_type}')

    # Some dataframes might have empty cells so let's exchange them with nans
    df = df.replace(r'^\s*$', np.nan, regex=True)
    if '2nd generation' in prot_baskets:
        annot_dict = read_basket_annotation_generation2(prot_baskets, basket_type)
    elif '3rd gen' in prot_baskets:
        annot_dict = read_basket_annotation_generation2(prot_baskets, basket_type)
    elif '4th gen' in prot_baskets:
        # annot_dict = basket_scoring.read_baskets_file_4th_gen(prot_baskets)
        annot_dict = read_basket_annotation_generation4(prot_baskets, data_type)
    else:
        annot_dict = read_basket_annotation_generation1(prot_baskets)

    if 'fp' in data_type:
        gene_df = df.index.to_frame()
    elif 'pp' in data_type:
        gene_df = df[['Gene names']].fillna("")

    df[[basket_type, f'{basket_type}_weights']] = gene_df.apply(map_identifier_list_to_baskets, annot_dict=annot_dict,
                                                                basket_type=basket_type,
                                                                with_weights=True,
                                                                axis=1, result_type="expand")

    return df, annot_dict


def map_identifier_list_to_baskets(identifier_list: pd.Series,
                                   annot_dict: Dict[str, str],
                                   basket_type: str,
                                   with_weights: bool) -> pd.Series:
    """
    Takes a list of semicolon separated identifiers and returns the baskets
    matching the first identifier with basket annotations and weights
    Input identifier list has to be pd.Series and if method used via apply it has to be of type dataframe
    """
    # TODO: throw error if no annot_dict given
    for identifier in identifier_list[0].split(';'):
        if with_weights:
            basket = pd.Series(["", ""])
        else:
            basket = ""
        if identifier in annot_dict:
            if with_weights:
                return pd.Series([annot_dict[identifier][basket_type], annot_dict[identifier]['weight']])
            else:
                return annot_dict[identifier][basket_type]
    return basket


def read_basket_annotation_generation1(prot_baskets: str) -> Dict[str, str]:
    """
    The annotation files consists of baskets as columns, with each row a gene in the basket.
    """
    basket_annotation = pd.read_csv(prot_baskets, sep='\t')

    for column in basket_annotation.columns:
        if basket_annotation[column].dtype != np.float64 and basket_annotation[column].dtype != np.int64:
            basket_annotation[column] = basket_annotation[column].str.strip()

    # transform to long format, where each row is a gene-basket pair
    basket_annotation = basket_annotation.melt(value_name='gene', var_name='basket').drop_duplicates()
    basket_annotation['weight'] = 1
    
    basket_annotation = basket_annotation[~basket_annotation['gene'].duplicated(keep='first')]
    basket_annotation = basket_annotation.set_index('gene')

    return create_identifier_to_basket_dict(basket_annotation)


def read_basket_annotation_generation2(prot_baskets: str,
                                       basket_type: str) -> Dict[str, str]:
    """
    The annotation files consists of rows where each row is a (gene,basket,weight) triplet.
    """
    basket_annotation = pd.read_excel(prot_baskets, sheet_name=basket_type)

    for column in basket_annotation.columns:
        if basket_annotation[column].dtype != np.float64 and basket_annotation[column].dtype != np.int64:
            basket_annotation[column] = basket_annotation[column].str.strip()

    # take only the first gene name if there is a group of genes
    basket_annotation['gene'] = basket_annotation['gene'].apply(lambda x: x.split(';')[0])
    
    basket_annotation = basket_annotation[~basket_annotation['gene'].duplicated(keep='first')]
    basket_annotation = basket_annotation.set_index('gene')
    return create_identifier_to_basket_dict(basket_annotation)


def create_identifier_to_basket_dict(basket_annotation: pd.DataFrame, data_type: Union[str, None] = 'fp',
                                     identifier_type: str = 'gene') -> Dict[str, str]:
    """
    collect all the baskets per gene in a dictionary of {'gene_name': 'basket1;basket2;...'}
    """
    if 'fp' in data_type:
        accepted_type = ['expression', 'kinase activity']
        # remove non fp types
        basket_annotation = basket_annotation[basket_annotation['LEVEL'].isin(accepted_type)]
        basket_annotation = basket_annotation.groupby([identifier_type]).agg(lambda x: ";".join(map(str, x)))
    elif 'pp' in data_type:
        accepted_type = ['phosphorylation', 'important phosphorylation', 'kinase activity']
        # remove non pp types
        basket_annotation = basket_annotation[basket_annotation['LEVEL'].isin(accepted_type)]
        basket_annotation = basket_annotation.groupby([identifier_type]).agg(lambda x: ";".join(map(str, x)))
    annot_dict = basket_annotation.to_dict('index')
    return annot_dict


def read_basket_annotation_generation4(prot_baskets: str, data_type: str):
    basket_annotation = pd.read_excel(prot_baskets)
    basket_annotation = utils.whitespace_remover(basket_annotation)
    basket_annotation['subbasket_level'] = basket_annotation['SUBBASKET'] + " - " + basket_annotation['LEVEL']  # basket_annotation['BASKET'] + " - " +
    basket_annotation['WEIGHT'] = basket_annotation['WEIGHT'].fillna(1)  # empty cell in WEIGHT column means weight = 1
    basket_annotation = basket_annotation.rename(
        {'BASKET': 'basket', 'SUBBASKET': 'sub_basket', 'WEIGHT': 'weight', 'GENE NAME': 'gene'}, axis=1)
    return create_identifier_to_basket_dict(basket_annotation, data_type)


"""
python3 -m bin.clinical_tools -c config_patients.json -i <input_tsv> -o <output_tsv>
"""
if __name__ == "__main__":
    import argparse
    import json
    import time

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True,
                        help="Absolute path to configuration file.")
    parser.add_argument("-i", "--input_file", required=True,
                        help="Absolute path to a tab separated file.")
    parser.add_argument("-o", "--output_file", required=True,
                        help="Absolute path to output file.")
    parser.add_argument("-t", "--data_type", default='fp',
                        help="Data type, either 'pp' or 'fp' (default: 'fp').")

    args = parser.parse_args()

    with open(args.config, "r") as f:
        configs = json.load(f)

    index_col = 'Gene names'
    if args.data_type == 'pp':
        index_col = 'Modified sequence'

    df = pd.read_csv(args.input_file, sep='\t', index_col=index_col)
    basket_file = configs["clinic_proc"]["prot_baskets"]

    # Start pipeline
    t0 = time.time()
    # TODO: fix.. outdated basket scheme
    for basket_type, baskets in zip(['basket', 'rtk'], [basket_file, basket_file]):
        df = prot_basket_annotation(df, baskets,
                                    data_type=args.data_type,
                                    basket_type=basket_type)

    df.to_csv(args.output_file, sep='\t')

    t1 = time.time()
    total = t1 - t0
    logger.info(f"Basket annotation finished in {total} seconds")
