from pathlib import Path

import bin.simsi as simsi
from bin.preprocess_tools import MQ_EVIDENCE_COLUMNS
from bin.data_loaders.simsi_tmt_loader import SimsiTMTLoader


def test_find_simsi_evidence_file():
    results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.05.15_CJ_maxlfq_minimal_test'
    simsi_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI'
    data_type = 'FP'
    simsi_evidence_file = simsi.find_simsi_evidence_file(results_folder, simsi_folder, data_type)
    assert simsi_evidence_file is not None


def test_load_simsi_evidence_file_archived():
    results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.05.15_CJ_maxlfq_minimal_test'
    simsi_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI'
    data_type = 'FP'
    evidence_files = ['/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches/Chordoma/Batch63_FP_Chordoma/combined/txt/evidence.txt',
                      '/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches/Chordoma/Batch64_FP_Chordoma/combined/txt/evidence.txt',
                      '/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches/Sarcoma/Batch1_FP_INFORM_MASTER/combined/txt/evidence.txt']
    loader = SimsiTMTLoader(evidence_files, results_folder, simsi_folder, data_type)
    all_batches = loader.load_data(MQ_EVIDENCE_COLUMNS)
    assert len(all_batches) == 3


def test_find_matching_summaries_folder():
    """
    Test that meta_input_file_FP.tsv in a regular summaries folder can be found
    """
    simsi_output_folder = Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI/simsi_output/FP')
    meta_input_file = Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.05.17_AhS_Batch106/meta_input_file_FP.tsv')
    summaries_folder = simsi.find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    assert summaries_folder == Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI/simsi_output/FP/summaries_2023.05.17_AhS_Batch106')


def test_find_matching_summaries_folder_archived():
    """
    Test that meta_input_file_FP.tsv inside a zip archive can be found
    """
    simsi_output_folder = Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI/simsi_output/FP')
    meta_input_file = Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.05.17_AhS_bis_Batch103/meta_input_file_FP.tsv')
    summaries_folder = simsi.find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    assert summaries_folder == Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI/simsi_output/FP/summaries_2023.05.11_CJ_batch103')


def test_find_matching_summaries_folder_not_found():
    """
    Test that trying to find an FP SIMSI search in the PP SIMSI folder will not return any results
    """
    simsi_output_folder = Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/SIMSI/simsi_output/PP')
    meta_input_file = Path('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.05.17_AhS_Batch106/meta_input_file_FP.tsv')
    summaries_folder = simsi.find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    assert summaries_folder is None