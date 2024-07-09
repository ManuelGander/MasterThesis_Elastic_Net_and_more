from pathlib import Path

import pytest

import bin.simsi as simsi


def test_get_correction_factor_files(correction_file_mapping_file, queue_file):
    assert simsi.get_correction_factor_files(['Batch2_PP_MASTER'], correction_file_mapping_file, queue_file, Path('.')) == ['UF291262_UE277617.txt']


def test_get_correction_factor_files_from_queue_file(correction_file_mapping_file, queue_file):
    assert simsi.get_correction_factor_files(['Batch3_PP_MASTER'], correction_file_mapping_file, queue_file, Path('.')) == ['UF291262_UE277617.txt']


def test_get_correction_factor_files_missing(correction_file_mapping_file, queue_file):
    with pytest.raises(ValueError, match=r".*Batch34_PP_MASTER.*"):
        simsi.get_correction_factor_files(['Batch34_PP_MASTER'], correction_file_mapping_file, queue_file, Path('.'))
    

@pytest.fixture
def queue_file(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    queue_file = tmp_path / "queue.csv"
    queue_file_content = """drug, MQ version, pre payload, post payload, threads, experiment, fasta file, raw folder, phospho, protease, mqpar, tmt corr factors, peptide fdr, protein fdr
Batch23_PP_INFORM,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,VA296455_UJ279751.txt,0.01,1
Batch22_PP_INFORM,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,VA296455_UJ279751.txt,0.01,1
Batch1_PP_INFORM_MASTER,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,UF291262_UE277617.txt,0.01,1
Batch2_PP_MASTER,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,UF291262_UE277617.txt,0.01,1
Batch3_PP_MASTER,1.6.12.0,cpAndGenerateMQpar,cleanPostSearchAndQCTopasPP,12,Patients,uniprot_proteome_up000005640_03112020.fasta,y/lumos_2/raw/,1,Trypsin/P,mqpar_base_1.6.12.0.xml,UF291262_UE277617.txt,0.01,1"""
    
    with open(queue_file, 'w') as f:
        f.write(queue_file_content)
    return queue_file


@pytest.fixture
def correction_file_mapping_file(tmp_path):
    tmp_path.mkdir(exist_ok=True)
    correction_file_mapping_file = tmp_path / "correction_factor_files.csv"
    correction_file_mapping_file_content = """experiment	correction_factor_file
Batch23_PP_INFORM	VA296455_UJ279751.txt
Batch22_PP_INFORM	VA296455_UJ279751.txt
Batch1_PP_INFORM_MASTER	UF291262_UE277617.txt
Batch2_PP_MASTER	UF291262_UE277617.txt"""
    
    with open(correction_file_mapping_file, 'w') as f:
        f.write(correction_file_mapping_file_content)
    return correction_file_mapping_file
