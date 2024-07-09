import pandas as pd
import numpy as np

import bin.TOPAS_scoring_functions as scoring_functions

def test_calculate_local_peptide_weight():
    # Sample input dataframe
    df = pd.DataFrame({
        'Modified sequence': ['A', 'B', 'C', 'D'],
        'pat_1': [10, 20, 30, 40],
        'pat_2': [5, 10, 15, np.nan],
        'Peptide count': [1, 2, 3, 4],
        'Site positions': ['X_S123', 'Y_T234', 'X_S123', 'Y_T234']
    })

    # Expected output dataframe
    expected_df = pd.DataFrame({
        'Modified sequence': ['A', 'B', 'C', 'D'],
        'pat_1': [10, 20, 30, 40],
        'pat_2': [5, 10, 15, np.nan],
        'Peptide count': [1, 2, 3, 4],
        'Site positions': ['X_S123', 'Y_T234', 'X_S123', 'Y_T234'],
        'weight_1': [1/4, 2/6, 3/4, 4/6],
        'weight_2': [1/4, 1, 3/4, np.nan],
    })

    # Call the function
    result_df = scoring_functions.calculate_psite_weights(df)

    # Assert that the result matches the expected output
    pd.testing.assert_frame_equal(result_df, expected_df, check_like=True)


def test_calculate_weighted_zscores():
    df = pd.DataFrame({
        'Modified sequence': ['A', 'B', 'C', 'D'],
        'pat_1': [10, 20, 30, 40],
        'pat_2': [5, 10, 15, np.nan],
        'Peptide count': [1, 2, 3, 4],
        'Site positions': ['X_S123', 'Y_T234', 'X_S123', 'Y_T234'],
        'weight_1': [1/4, 2/6, 3/4, 4/6],
        'weight_2': [1/4, 1, 3/4, np.nan],
    })

    expected_df = pd.DataFrame({
        'Modified sequence': ['A', 'B', 'C', 'D'],
        'pat_1': [10, 20, 30, 40],
        'pat_2': [5, 10, 15, np.nan],
        'Peptide count': [1, 2, 3, 4],
        'Site positions': ['X_S123', 'Y_T234', 'X_S123', 'Y_T234'],
        'weight_1': [1/4, 2/6, 3/4, 4/6],
        'weight_2': [1/4, 1, 3/4, np.nan],
        'weighted_1': [2.5, 20/3, 22.5, 80/3],
        'weighted_2': [5/4, 10, 45/4, np.nan],
    })

    # Call the function
    result_df = scoring_functions.calculate_weighted_z_scores(df)

    # Assert that the result matches the expected output
    pd.testing.assert_frame_equal(result_df, expected_df, check_like=True)