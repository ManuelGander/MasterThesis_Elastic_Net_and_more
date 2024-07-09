import re
import toml
import os

import pandas as pd
import matplotlib.pyplot as plt

from ... import drug_scoring

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

folder = '2022.08.02_MT_mixed_cohort'
#folder = '2022.08.02_MT_mixed_cohort_with_SIMSI'

wp3_df = pd.read_csv(f'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/{folder}/annot_pp.csv')
wp3_df['Modified sequence'] = wp3_df['Modified sequence'].apply(update_modifications)
wp3_df['Sequence'] = wp3_df['Modified sequence'].apply(remove_modifications)

wp1_df = drug_scoring.load_topas_wp1_data('all')
wp1_df['Sequence'] = wp1_df['Modified sequence'].apply(remove_modifications)


print("Modified sequences WP1:", len(set(wp1_df['Modified sequence'])))
print("Modified sequences WP3:", len(set(wp3_df['Modified sequence'])))
print("Modified sequences overlap:", len(set(wp3_df['Modified sequence']).intersection(set(wp1_df['Modified sequence']))))

print("Sequences WP1:", len(set(wp1_df['Sequence'])))
print("Sequences WP3:", len(set(wp3_df['Sequence'])))
print("Sequences overlap:", len(set(wp3_df['Sequence']).intersection(set(wp1_df['Sequence']))))

def update_modifications(mod_sequence):
    raw_sequence = mod_sequence.replace('pS', 'S(ph)')
    raw_sequence = raw_sequence.replace('pT', 'T(ph)')
    raw_sequence = raw_sequence.replace('pY', 'Y(ph)')
    raw_sequence = raw_sequence.replace('(Acetyl (Protein N-term))', '')
    raw_sequence = raw_sequence.replace('(Oxidation (M))', '')
    raw_sequence = raw_sequence.replace('_', '')
    return raw_sequence


def remove_modifications(mod_sequence):
    raw_sequence = mod_sequence.replace('S(ph)', 'S')
    raw_sequence = raw_sequence.replace('T(ph)', 'T')
    raw_sequence = raw_sequence.replace('Y(ph)', 'Y')
    return raw_sequence
