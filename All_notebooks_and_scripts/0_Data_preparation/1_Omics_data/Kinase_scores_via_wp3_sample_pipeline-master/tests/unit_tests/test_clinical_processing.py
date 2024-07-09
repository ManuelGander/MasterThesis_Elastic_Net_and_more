import pytest
import pandas as pd

import bin.clinical_tools as ct


def test_read_basket_annotation_generation1():
    prot_baskets = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches/Sarcoma/basketannotation_kegg_wiki_MASTER.txt"
    basket_type = "N/A"
    
    annot_dict = ct.read_basket_annotation_generation1(prot_baskets)
    print(annot_dict)
    assert annot_dict['ZNRF3']['basket'] == 'Wnt Notch Hedgehog'


def test_read_basket_annotation_generation2():
    prot_baskets = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Baskets_RTK Baskets_2nd generation_220601.xlsx"
    basket_type = "basket"
    
    annot_dict = ct.read_basket_annotation_generation2(prot_baskets, basket_type)
    assert annot_dict['WEE1']['basket'] == 'Cell Cycle'
    assert annot_dict['WEE1']['weight'] == -1


def test_prot_basket_annotation_generation1():
    prot_baskets = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches/Sarcoma/basketannotation_kegg_wiki_MASTER.txt"
    basket_type = "basket"
    data_type = "pp"
    df = pd.DataFrame({'Gene names': ['ZNRF3']})
    
    basket_annotated_df, _ = ct.prot_basket_annotation(df, prot_baskets, data_type, basket_type)    
    assert basket_annotated_df.loc[0]['basket'] == 'Wnt Notch Hedgehog'


def test_prot_basket_annotation_generation2():
    prot_baskets = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Baskets_RTK Baskets_2nd generation_220601.xlsx"
    basket_type = "basket"
    data_type = "pp"
    df = pd.DataFrame({'Gene names': ['WEE1']})
    
    basket_annotated_df, _ = ct.prot_basket_annotation(df, prot_baskets, data_type, basket_type)    
    assert basket_annotated_df.loc[0]['basket'] == 'Cell Cycle'


def test_prot_basket_annotation_generation2_fp():
    prot_baskets = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Baskets_RTK Baskets_2nd generation_220601.xlsx"
    basket_type = "basket"
    data_type = "fp"
    df = pd.DataFrame({'Gene names': ['WEE1']}).set_index('Gene names')
    
    basket_annotated_df, _ = ct.prot_basket_annotation(df, prot_baskets, data_type, basket_type)
    assert basket_annotated_df.loc['WEE1']['basket'] == 'Cell Cycle'

