import io

import pytest
import numpy as np
import pandas as pd

from bin.preprocess_tools import impute_data, convert_long_to_wide_format, filter_evidence_files, median_centering
from bin.sample_annotation import get_unique_batches
from bin.identification_metadata import mark_num_peptides
from bin.data_loaders import data_loader, tmt_loader

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


class TestImputation:

    """
    we need minimum 1 channel for imputation
    # check that minimum is chosen when minimum and max/100 when that's the case:   DONE
    # check imputation status added:                                                DONE
    # check that only patient channels are imputed:                                 DONE
    # check that there are no 0s:
    # perhaps just some rows for each test making it clear what is tested and what breaks !

    """

    # TODO: use fixture with readable dataframe
    def test_impute_data(self, df_tmt11, exp_df_tmt11_imputed_with_status):
        test_df = impute_data(df_tmt11)
        exp_result = exp_df_tmt11_imputed_with_status
        # hacky way to test only imputation values
        test_df = test_df.loc[:, test_df.columns.str.startswith('Reporter')]
        exp_result = exp_result.loc[:, exp_result.columns.str.startswith('Reporter')]
        pd.testing.assert_frame_equal(test_df, exp_result, atol=.01)

    def test_impute_data_imputation_status(self, df_tmt11, exp_df_tmt11_imputed_with_status):
        test_df = impute_data(df_tmt11)
        exp_result = exp_df_tmt11_imputed_with_status

        # hacky way to test only imputation values
        test_df = test_df.loc[:, (test_df.columns.str.startswith('Reporter') | test_df.columns.str.startswith('Identification'))]
        pd.testing.assert_frame_equal(test_df, exp_result, atol=.01)


def test_impute_data_missing_ref(self):
    test_df = pd.DataFrame(
        [(19022.0, 14528.0, 32250.0, 23406.0, 621154.0, 9854.9, 141170.00, 71932.0, 90506.0, 57724.0, 57168.0,
          'x', '', '', '', '', '', '', '', '', '', ''),
         (np.nan, 343490.0, np.nan, np.nan, 310470.0, 336740.0, 550.37, np.nan, 72777.0, 57348.0, 51079.0,
          'x;', '', '', '', '', '', '', '', '', '', ''),
         (50365.0, 18493.0, 69062.0, np.nan, 16222.0, 5868.5, 3980.80, 5669.2, np.nan, np.nan, np.nan,
          '', '', '', '', '', '', '', '', '', '', ''),
         (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
          '', '', '', '', '', '', '', '', '', '', '')],
        columns=[f'Reporter intensity corrected {i}' for i in range(1, 12)] + [f'Identification metadata {i}' for i in range(1, 12)])
    # print(test_df)
    test_df = impute_data(test_df)
    expect_result = pd.DataFrame(
        [(19022.0, 14528.0, 32250.0, 23406.0, 621154.0, 9854.9, 141170.00, 71932.0, 90506.0, 57724.0, 57168.0,
          'x', '', '', '', '', '', '', '', '', '', ''),
         (550.37, 343490.0, 550.37, 550.37, 310470.0, 336740.0, 550.37, 550.37, 72777.0, 57348.0, 51079.0,
          'x;imputed;', '', 'imputed;', 'imputed;', '', '', '', 'imputed;', '', '', ''),
         (50365.0, 18493.0, 69062.0, 690.62, 16222.0, 5868.5, 3980.80, 5669.2, np.nan, np.nan, np.nan,
          '', '', '', 'imputed;', '', '', '', '', '', '', ''),
         (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
          '', '', '', '', '', '', '', '', '', '', '')],
        columns=[f'Reporter intensity corrected {i}' for i in range(1, 12)] + [f'Identification metadata {i}' for i in range(1, 12)])
    pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)


class TestMedianCentering:
    # TODO: make simpler input dataframe where we can easily see the median
    def test_median_centering(self):
        test_df = pd.DataFrame(
            [(1, 1, 1, 1, 1, 1, 1, 1, 4, 3, 2),
             (1, 2, 1, 2, 1, 1, 1, 1, 4, 3, 2),
             (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
             (2, 2, 2, 3, 6, 2, 2, 2, 5, 6, 3)])
        # average of the column medians is 2.23
        test_df = median_centering(test_df)
        expect_result = pd.DataFrame(
            [[1.48, 1.11, 1.48, 1.11, 1.48, 1.48, 1.48, 1.48, 1.98, 1.48, 1.78],
             [1.48, 2.23, 1.48, 2.23, 1.48, 1.48, 1.48, 1.48, 1.98, 1.48, 1.78],
             [2.97, 2.23, 2.97, 2.23, 2.97, 2.97, 2.97, 2.97, 2.47, 2.97, 2.67],
             [2.97, 2.23, 2.97, 3.34, 8.91, 2.97, 2.97, 2.97, 2.47, 2.97, 2.67]])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)

    def test_median_centering_with_zeroes(self):
        test_df = pd.DataFrame(
            [(0, 1, 1, 1, 6, 1, 1, 1, 4, 3, 2),
             (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
             (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2),
             (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3)])
        # average of the column medians is 2.318
        test_df = median_centering(test_df)
        expect_result = pd.DataFrame(
            [[np.nan, 1.16, 1.55, 1.16, 6.95, 1.55, 1.55, 1.55, 2.06, 1.55, 1.85],
             [2.32, 2.32, 3.09, 2.32, 2.32, 3.09, 3.09, 3.09, 2.58, 3.09, 2.78],
             [1.16, 2.32, 1.55, 3.48, 1.16, 1.55, 1.55, 1.55, 2.06, 1.55, 1.85],
             [2.32, 2.32, 3.09, 2.32, 2.32, 3.09, 3.09, 3.09, 2.58, 3.09, 2.78]])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)

    def test_median_centering_too_many_zeroes_in_row(self):
        test_df = pd.DataFrame(
            [(0, 0, 0, 0, 6, 1, 1, 1, 4, 3, 2),
             (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
             (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2),
             (2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3)])
        # average of the column medians is 2.318
        test_df = median_centering(test_df)
        expect_result = pd.DataFrame(
            [[np.nan, np.nan, np.nan, np.nan, 8.18, 1.36, 1.36, 1.36, 2.18, 1.36, 1.82],
             [2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73],
             [1.36, 2.73, 1.36, 4.09, 1.36, 1.36, 1.36, 1.36, 2.18, 1.36, 1.82],
             [2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73, 2.73]])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)


#    def test_sample_correction(self):
#        # check that 0 turns into nan
#        test_df = pd.DataFrame(
#            [(1, 1, 1, 1, 6, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
#             (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3)])
#        test_df = preprocess.sample_correction(test_df)
#        print(test_df)
#        expect_result = pd.DataFrame(
#            [(1, 1, 1, 1, 6, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3),
#             (1, 2, 1, 3, 1, 1, 1, 1, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 6, 3)])
#        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)


class TestRedistributeMS1:
    def test_ms1_imputation(self, long_df_tmt11):
        df = tmt_loader._impute_ms1_intensity(long_df_tmt11)
        # print(df)
        assert df.iloc[1]['MS1'] == 5.0  # missing MS1, same channel intensities
        assert df.iloc[2]['MS1'] == 5.0  # missing MS1, doubled all channel intensities ==> should not have any influence
        assert df.iloc[5]['MS1'] == 6.0  # missing MS1, average of other MS1s if QC channels have the same ratio relative to the summed
        # intensity
        assert df.iloc[6]['MS1'] == 3.0  # missing MS1, doubled intensity in reference channel, sum of all channels unchanged
        assert np.isnan(df.iloc[8]['MS1'])  # missing MS1, no other batches to transfer from

    def test_redistribute_ms1_intensity(self, long_df_tmt11):
        df = tmt_loader._impute_ms1_intensity(long_df_tmt11)
        df = tmt_loader._redistribute_ms1_intensity(df)
        # print(df)
        assert df.iloc[0]['Reporter intensity corrected 1'] == 0.25  # 1.0 (TMT) / 20.0 (TMT sum) * 5.0 (MS1) = 0.25
        assert df.iloc[1]['Reporter intensity corrected 1'] == 0.25  # 1.0 / 20.0 * 5.0 = 0.25
        assert df.iloc[2]['Reporter intensity corrected 1'] == 0.25  # 1.0 / 20.0 * 5.0 = 0.25
        assert np.isnan(df.iloc[5]['Reporter intensity corrected 2'])  # do not change abundance if TMT abundance is missing
        assert df.iloc[8]['Reporter intensity corrected 1'] == 3.00  # do not change abundance if MS1 is missing


def test_convert_long_to_wide_format(long_df):
    # print(long_df)
    wide_df = convert_long_to_wide_format(long_df)
    # print(wide_df)
    assert len(wide_df.index) == 4
    assert len(wide_df.dropna().index) == 2
    # TODO: add more assertions to check for gene name and protein aggregation


class TestIdentificationMetadata:
    """

    Check that there can already be imputed values in there before adding num_peptides: DONE
    Check that 0 does not give number of peptides: DONE
    Check that nan does not give number of peptides:
    Check that we get ; after:
    TODO: use fixture. add row with nan example

    Do we want to fix the xnum_peptides situations?
    """

    def test_mark_num_peptides(self):
        test_df = pd.DataFrame(
            [(2, 0, 2, 2, 2, 2, 1, 2, 2, 2, 3,
              'x', '', '', '', '', '', '', '', '', '', ''),
             (3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
              'x;', '', '', '', '', '', '', '', '', '', 'y;')],
            columns=[f'Unique peptides Reporter intensity corrected {i}' for i in range(1, 12)] + [f'Identification metadata {i}' for i in
                                                                                                   range(1, 12)])
        test_df = mark_num_peptides(test_df)
        print(test_df)
        expect_result = pd.DataFrame(
            [(2, 0, 2, 2, 2, 2, 1, 2, 2, 2, 3,
              'xnum_peptides=2;', '', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;',
              'num_peptides=1;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=3;'),
             (3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
              'x;num_peptides=3;', 'num_peptides=1;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;',
              'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'y;num_peptides=1;')],
            columns=[f'Unique peptides Reporter intensity corrected {i}' for i in range(1, 12)] + [f'Identification metadata {i}' for i in
                                                                                                   range(1, 12)])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)

    def test_mark_num_peptides_with_na(self):
        test_df = pd.DataFrame(
            [(2, np.nan, 2, 2, 2, 2, 1, 2, 2, 2, 3,
              'x', '', '', '', '', '', '', '', '', '', ''),
             (3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
              'x;', '', '', '', '', '', '', '', '', '', 'y;')],
            columns=[f'Unique peptides Reporter intensity corrected {i}' for i in range(1, 12)] + [f'Identification metadata {i}' for i in
                                                                                                   range(1, 12)])
        test_df = mark_num_peptides(test_df)
        print(test_df)
        expect_result = pd.DataFrame(
            [(2, np.nan, 2, 2, 2, 2, 1, 2, 2, 2, 3,
              'xnum_peptides=2;', '', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;',
              'num_peptides=1;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=3;'),
             (3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1,
              'x;num_peptides=3;', 'num_peptides=1;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;',
              'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'num_peptides=2;', 'y;num_peptides=1;')],
            columns=[f'Unique peptides Reporter intensity corrected {i}' for i in range(1, 12)] + [f'Identification metadata {i}' for i in
                                                                                                   range(1, 12)])
        pd.testing.assert_frame_equal(test_df, expect_result, atol=.01)


class TestBatchHandling:
    """

    """

    def test_get_unique_batches(self):
        batches = pd.DataFrame(
            [(1, 'Sarcoma'),
             (1, 'Sarcoma'),
             (1, 'Sarcoma'),
             (6, 'Sarcoma'),
             (1, 'Chordoma'),
             (1, 'Chordoma'),
             (1, 'Chordoma'),
             (4, 'Chordoma')])
        batches.columns = ['Batch Name', 'Cohort']
        assert get_unique_batches(batches) == [[1, 'Sarcoma'], [6, 'Sarcoma'], [1, 'Chordoma'], [4, 'Chordoma']]

    def test_filter_evidence_files(self):
        evidence_files = ['/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch2_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch3_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch4_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch5_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Chordoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Chordoma/Batch4_FP_Blabla/combined/txt/evidence.txt']
        batches = [[1, 'Sarcoma'], [6, 'Sarcoma'], [1, 'Chordoma'], [4, 'Chordoma']]
        filtered_evidence_files = filter_evidence_files(evidence_files, 'FP', batches)
        assert filtered_evidence_files == ['/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                                           '/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt',
                                           '/my/path/Chordoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                                           '/my/path/Chordoma/Batch4_FP_Blabla/combined/txt/evidence.txt']

    def test_extract_batch_name(self):
        evidence_files = ['/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Chordoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Chordoma/Batch14_FP_Blabla/combined/txt/evidence.txt']

        extracted_batches = [data_loader.extract_batch_name(e) for e in evidence_files]
        assert extracted_batches == ['Sarcoma_Batch1', 'Sarcoma_Batch6', 'Chordoma_Batch1', 'Chordoma_Batch14']

    def test_extract_experiment_name(self):
        evidence_files = ['/my/path/Sarcoma/Batch1_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Sarcoma/Batch6_FP_Blabla/combined/txt/evidence.txt',
                          '/my/path/Chordoma/Batch1_FP_Blabla2/combined/txt/evidence.txt',
                          '/my/path/Chordoma/Batch14_FP_Blabla/combined/txt/evidence.txt']

        extracted_batches = [data_loader.extract_experiment_name(e) for e in evidence_files]
        assert extracted_batches == ['Batch1_FP_Blabla', 'Batch6_FP_Blabla', 'Batch1_FP_Blabla2', 'Batch14_FP_Blabla']


# Creating dataframes from strings: https://towardsdatascience.com/67b0c2b71e6a
@pytest.fixture
def df_tmt11():
    df_string = """Sequence,                         Modified sequence,    Gene names,      Proteins, Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3, Reporter intensity corrected 4, Reporter intensity corrected 5, Reporter intensity corrected 6, Reporter intensity corrected 7, Reporter intensity corrected 8, Reporter intensity corrected 9, Reporter intensity corrected 10, Reporter intensity corrected 11
DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_A;GENE_B, PROT_A;PROT_B, , , , , , , , , , , 
DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0
DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A, , , , , , , , , 6.0, 2.0, 2.0
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_C;GENE_B, PROT_A;PROT_B, 3.0, , , , , , , , , ,
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 3.0, 3.0,    , 3.0, 3.0,    , 1.0, 2.0,    , 3.0, 2.0
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 3.0,    , 3.0, 3.0,    , 3.0, 1.0, 2.0,    , 3.0, 2.0
     SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 0.03,    , 3.0, 3.0,    , 3.0, 1.0, 1.0,    , 6.0,    
     QSPTPESPTMLTK,      QS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, , , , , , , , , 1.0, 3.0, 2.0
     KSPTPESPTMLTK,      KS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, 3.0, 2.0, , , , 1.0, 3.0, 2.0, , , 
"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=',', skipinitialspace=True)
    for i in range(1, 12):
        df[f'Identification metadata {i}'] = ""
    return df


@pytest.fixture
def exp_df_tmt11_imputed_with_status():
    df_string = """ Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3, Reporter intensity corrected 4, Reporter intensity corrected 5, Reporter intensity corrected 6, Reporter intensity corrected 7, Reporter intensity corrected 8, Reporter intensity corrected 9, Reporter intensity corrected 10, Reporter intensity corrected 11, Identification metadata 1, Identification metadata 2, Identification metadata 3, Identification metadata 4, Identification metadata 5, Identification metadata 6, Identification metadata 7, Identification metadata 8, Identification metadata 9, Identification metadata 10, Identification metadata 11
NaN,  NaN,   NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN, NaN, '', '', '', '', '', '', '', '', '', '', ''
1.00, 2.00, 3.00, 1.00, 2.00, 3.00, 1.00, 2.00,  3.0,  1.0, 1.0, '', '', '', '', '', '', '', '', '', '', ''
0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,  6.0,  2.0, 2.0, 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', '', '', ''
3.00, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,  NaN,  NaN, NaN, '', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', '', '', ''
3.00, 3.00, 0.03, 3.00, 3.00, 0.03, 1.00, 2.00,  NaN,  3.0, 2.0, '', '', 'imputed;', '', '', 'imputed;', '', '', '', '', ''
3.00, 0.03, 3.00, 3.00, 0.03, 3.00, 1.00, 2.00,  NaN,  3.0, 2.0, '', 'imputed;', '', '', 'imputed;', '', '', '', '', '', ''
0.03, 0.03, 3.00, 3.00, 0.03, 3.00, 1.00, 1.00,  NaN,  6.0, NaN, '', 'imputed;', '', '', 'imputed;', '', '', '', '', '', ''
0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,  1.0,  3.0, 2.0, 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', 'imputed;', '', '', ''
3.00, 2.00, 0.03, 0.03, 0.03, 1.00, 3.00, 2.00,  NaN,  NaN, NaN, '', '', 'imputed;', 'imputed;', 'imputed;', '', '', '', '', '', ''
    """
    df = pd.read_csv(io.StringIO(df_string),  dtype={f'Reporter intensity corrected {i}': float for i in range(1, 12)}, delimiter=',', skipinitialspace=True, quotechar="'",
                     keep_default_na=False, na_values='NaN')
    return df


@pytest.fixture
def long_df():
    df_string = """Batch, Sequence,                         Modified sequence,    Gene names,      Proteins, Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3
Batch_1, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_A;GENE_B, PROT_A;PROT_B, 0.5, 5.3, 9.2
Batch_1,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_C;GENE_B, PROT_A;PROT_B, 0.4, 7.3, 9.3
Batch_1, WSDSWDADAFSVEDPVRK, WS(Phospho (STY))DSWDADAFS(Phospho (STY))VEDPVRK, GENE_D;GENE_B, PROT_A;PROT_B, 0.3, 1.3, 9.4
Batch_2, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A, 0.6, 2.3, 9.6
Batch_2,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 0.8, 1.3, 9.8
Batch_2,      TSPTPESPTMLTK,      TS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, 0.7, 9.3, 9.7"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=',', skipinitialspace=True)
    df['Modifications'] = ""
    df['Protein Names'] = ""
    df['Charge'] = 2
    df['m/z'] = 500.0
    df['Mass'] = 1000.0
    df['Missed cleavages'] = 0
    df['Length'] = 1
    df['Reverse'] = ""
    for i in range(1, 4):
        df[f'Identification metadata {i}'] = ""
    return df


@pytest.fixture
def long_df_tmt11():
    df_string = """Batch, Sequence,                         Modified sequence,    Gene names,      Proteins, MS1, Reporter intensity corrected 1, Reporter intensity corrected 2, Reporter intensity corrected 3, Reporter intensity corrected 4, Reporter intensity corrected 5, Reporter intensity corrected 6, Reporter intensity corrected 7, Reporter intensity corrected 8, Reporter intensity corrected 9, Reporter intensity corrected 10, Reporter intensity corrected 11
Batch_1, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_A;GENE_B, PROT_A;PROT_B, 5.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0
Batch_2, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A,    , 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0
Batch_2, DSDSWDADAFSVEDPVRK, DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK, GENE_B;GENE_A, PROT_B;PROT_A,    , 2.0, 4.0, 6.0, 2.0, 4.0, 6.0, 2.0, 4.0, 6.0, 2.0, 2.0
Batch_1,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_C;GENE_B, PROT_A;PROT_B, 4.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 1.0, 1.0, 2.0, 4.0,    
Batch_2,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B, 8.0, 3.0, 3.0,    , 3.0, 3.0,    , 1.0, 2.0,    , 3.0, 2.0
Batch_3,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B,    , 3.0,    , 3.0, 3.0,    , 3.0, 1.0, 2.0,    , 3.0, 2.0
Batch_4,      SSPTPESPTMLTK,      SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK,        GENE_B,        PROT_B,    , 3.0,    , 3.0, 3.0,    , 3.0, 1.0, 1.0,    , 6.0,    
Batch_2,      TSPTPESPTMLTK,      TS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B, 0.7, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0
Batch_3,      QSPTPESPTMLTK,      QS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B,    , 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0
Batch_3,      KSPTPESPTMLTK,      KS(Phospho (STY))PTPES(Phospho (STY))PTMLTK, GENE_A;GENE_B,        PROT_B,    , 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0,    ,    , 2.0
"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=',', skipinitialspace=True)
    df['Modifications'] = ""
    df['Protein Names'] = ""
    df['Charge'] = 2
    df['m/z'] = 500.0
    df['Mass'] = 1000.0
    df['Missed cleavages'] = 0
    df['Length'] = 1
    df['Reverse'] = ""
    for i in range(1, 4):
        df[f'Identification metadata {i}'] = ""
    return df
