import io
import os
import sys
import argparse
from pathlib import Path
import subprocess
import logging
import time
import json
from typing import List, Dict, Union
import traceback
import zipfile

import pandas as pd
import simsi_transfer.main as simsi
from job_pool import JobPool

import bin.preprocess_tools as prep
from . import config
from . import sample_annotation
from .utils import init_file_logger, send_slack_message

# hacky way to get the package logger instead of just __main__ when running as python -m bin.simsi ...
logger = logging.getLogger(__package__ + "." + __file__)


def main(argv):

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", dest="config", required=True,
                        help="Absolute path to configuration file.")
    args = parser.parse_args(argv)

    configs = config.load(args.config)

    if not configs["preprocessing"]["run_simsi"]:
        logger.info(f"run_simsi flag is set to False, skipping SIMSI")
        return

    # Check if we have already processed results in our results folder and skip if so
    # data_types = configs["data_types"]
    # for data_type in data_types:
    #     flag = False
    #     if os.path.exists(os.path.join(configs["results_folder"], f'preprocessed_{data_type}2.csv')):
    #         flag = True
    #     if not flag:
    #         break
    # if flag:
    #     logger.info(f"Found already processed results in results folder, skipping SIMSI")
    #     return

    os.makedirs(configs["results_folder"], exist_ok=True)
    
    jsonString = json.dumps(configs, indent=4)
    with open(f'{configs["results_folder"]}/configs.json', "w") as jsonFile:
        jsonFile.write(jsonString)

    # Check for sample + metadata annotation errors before starting pipeline
    _ = prep.check_annot(configs["sample_annotation"], configs["metadata_annotation"], prep.in_metadata, configs["patient_regex"])

    # check_config(configs)

    run_simsi(configs["results_folder"], configs["preprocessing"]["raw_data_location"], configs["sample_annotation"],
              **configs["simsi"], data_types=configs["data_types"])


# def check_config(configs):
#
#     # check data types
#     for data_type in configs["data_types"]:
#         if data_type.upper() not in ['FP', 'PP']:
#             raise ValueError(f'Data type {data_type} not accepted. Accepted data types are `fp` (full proteome) and `pp` (phospho).')
#
#     # check existence of file locations
#     if not os.path.exists(configs["sample_annotation"]):
#         raise ValueError(f'Sample annotation file {configs["sample_annotation"]} does not exist in this location.')
#     metadata_file = configs["metadata_annotation"]
#     if not os.path.isfile(metadata_file):
#             metadata_file = metadata_file.replace('Retrospective_MTBs_Evaluation',
#                                                   'Retrospective_MTBs_Evaluation/Metadata_excel_outdated')
        # if not os.path.isfile(metadata_file):
        #     raise ValueError(f'Metadata: {metadata_file} not found.')





def run_simsi(*args, **kwargs) -> None:
    # propagate logs from simsi_transfer package to our current log handlers
    logging.getLogger('simsi_transfer').handlers = logging.getLogger(__package__).handlers
    
    data_types = kwargs.pop('data_types')
    processingPool = JobPool(processes=2, timeout=60000)  # 60,000 seconds = 16.7 hours
    for data_type in data_types:
        kwargs_with_data_type = kwargs.copy()
        kwargs_with_data_type['data_type'] = data_type

        processingPool.applyAsync(run_simsi_data_type, args, kwargs_with_data_type)
    processingPool.checkPool(printProgressEvery=1)


def run_simsi_data_type(results_folder: str,
              raw_data_location: str,
              sample_annotation_file: str,
              simsi_folder: str,
              tmt_ms_level: str,
              stringencies: int,
              tmt_requantify: bool,
              maximum_pep: int,
              num_threads: int,
              data_type: str):
    init_file_logger(results_folder, f'SIMSI_log_{data_type}.txt')
        
    # Start pipeline
    t0 = time.time()
    
    # run simsi (~10 hours)    
    logger.info(f"SIMSI started")
    
    try:
        summary_files = get_summary_files(raw_data_location, data_type, sample_annotation_file)
        result_folder_name = Path(results_folder).name
        
        copy_files(summary_files, data_type, simsi_folder)

        meta_input_file = get_meta_input_file_path(results_folder, data_type)
        create_meta_input_file(summary_files, data_type, simsi_folder, meta_input_file)

        simsi_output_folder = get_simsi_output_folder(simsi_folder, data_type)
        if matched_summaries_folder := find_matching_summaries_folder(simsi_output_folder, meta_input_file):
            logger.info(
                f"Found a summaries folder that matches the current list of folders: {matched_summaries_folder}\nSkipping SIMSI processing.")
            return

        run_simsi_single(meta_input_file, simsi_output_folder, tmt_ms_level, stringencies, tmt_requantify, maximum_pep, num_threads)

        store_results_for_reuse(simsi_output_folder, meta_input_file, result_folder_name)
        # empty_raw_files(simsi_folder, meta_input_file)
        message = f'SIMSI {data_type} finished'
    except Exception as e:
        logger.info(str(e))
        logger.info(traceback.format_exc())
        message = str(e)

    send_slack_message(message)

    t1 = time.time()
    total = t1 - t0
    logger.info(f"SIMSI finished in {total} seconds")


# def empty_raw_files(simsi_folder, meta_input_file):
#     # subtract one layer and go to raw files
#     simsi_raw_files_folder = Path(simsi_folder) / Path('raw_files')

#     # find list of what to empty

#     # only empty if they are not already?


def get_simsi_output_folder(simsi_folder: str, data_type: str):
    return Path(simsi_folder) / Path('simsi_output') / Path(data_type.upper())


def get_meta_input_file_path(results_folder: str, data_type: str):
    return Path(results_folder) / Path(f'meta_input_file_{data_type.upper()}.tsv')


def find_matching_summaries_folder(simsi_output_folder: Path, meta_input_file: Path):
    """
    Iterates over folders starting with 'summaries_' to see if any of them contains results for the current list of folders
    """
    summary_folders_and_archives = list(simsi_output_folder.glob('summaries_*'))
    summary_folders_and_archives.sort(key=os.path.getmtime, reverse=True)
    for s in summary_folders_and_archives:
        if s.is_dir():
            meta_input_file_other = s / Path(meta_input_file.name)
            if not (meta_input_file.is_file() and meta_input_file_other.is_file()):
                continue
        elif s.name.endswith(".zip"):
            try:
                archive = zipfile.ZipFile(s, 'r')
            except zipfile.BadZipFile:
                logger.warning(f"Encountered corrupt zip file: {s}")

            file_path_in_archive = s.with_suffix('').name + '/' + meta_input_file.name
            if not file_path_in_archive in archive.namelist():
                continue
            meta_input_file_other = io.BytesIO(archive.read(file_path_in_archive))
            s = s.with_suffix('')
        else:
            continue
        
        if meta_input_files_equal(meta_input_file, meta_input_file_other):
            return s

    return None


def meta_input_files_equal(meta_input_file: Path, meta_input_file_other: Path):
    meta_input_df = pd.read_csv(meta_input_file, sep='\t', usecols=['mq_txt_folder', 'raw_folder'])
    meta_input_other_df = pd.read_csv(meta_input_file_other, sep='\t', usecols=['mq_txt_folder', 'raw_folder'])
    return meta_input_df.sort_values(by='mq_txt_folder').equals(meta_input_other_df.sort_values(by='mq_txt_folder'))


def get_correction_factor_files(experiments: List[str], correction_file_mapping_file: Path, mq_queue_file: Path,
                                correction_file_folder: Path):
    # has two columns: experiment, correction_factor_file
    correction_file_mapping_df = pd.read_csv(str(correction_file_mapping_file), sep='\t', index_col='experiment')
    correction_files = map_experiments_to_files(experiments, correction_file_mapping_df, 'correction_factor_file')
    if None in correction_files:
        missing_experiments = get_missing_experiments(experiments, correction_files)
        update_correction_factor_file_mapping(missing_experiments, correction_file_mapping_file, mq_queue_file)
        # TODO: following caused error due to missing last input (now added for test?)
        return get_correction_factor_files(experiments, correction_file_mapping_file, mq_queue_file, correction_file_folder)

    return [str(correction_file_folder / f) for f in correction_files]


def update_correction_factor_file_mapping(missing_experiments: List[str], correction_file_mapping_file: Path,
                                          mq_queue_file: Path) -> None:
    # has many columns, including: drug (=experiment), ' tmt corr factors' (=correction_factor_file)
    queue_df = pd.read_csv(str(mq_queue_file), index_col='drug')
    correction_files_new = map_experiments_to_files(missing_experiments, queue_df, ' tmt corr factors')

    if None in correction_files_new:
        missing_experiments_new = get_missing_experiments(missing_experiments, correction_files_new)
        raise ValueError(f"Could not find experiments {missing_experiments_new} in Queue file")

    correction_file_mapping_df = pd.read_csv(str(correction_file_mapping_file), sep='\t')
    correction_file_mapping_df_new = pd.DataFrame(
        {'experiment': missing_experiments, 'correction_factor_file': correction_files_new})
    correction_file_mapping_df = pd.concat([correction_file_mapping_df, correction_file_mapping_df_new])

    correction_file_mapping_df.to_csv(str(correction_file_mapping_file), sep='\t', index=False)


def get_missing_experiments(experiments: List[str], correction_files: List[str]):
    return [experiment for experiment, correction_file in zip(experiments, correction_files) if correction_file is None]


def map_experiments_to_files(experiments: List[str], experiment_mapping_df: pd.DataFrame, key: str):
    experiment_mapping_df = drop_duplicate_indices(experiment_mapping_df)
    experiment_mapping = experiment_mapping_df.to_dict('index')
    return [experiment_mapping.get(experiment, {key: None})[key] for experiment in experiments]


def drop_duplicate_indices(df):
    return df[~df.index.duplicated(keep='last')]


def find_simsi_evidence_file(results_folder: str, simsi_folder: str, data_type: str) -> Union[str, io.BytesIO]:
    meta_input_file = get_meta_input_file_path(results_folder, data_type.upper())
    simsi_output_folder = get_simsi_output_folder(simsi_folder, data_type.upper())
    simsi_summaries_folder = find_matching_summaries_folder(simsi_output_folder, meta_input_file)
    if not simsi_summaries_folder:
        raise ValueError("Could not find a SIMSI summaries folder with a matching list of input folders")

    simsi_evidence_file = f'{str(simsi_summaries_folder)}/p10/p10_evidence.txt'
    logger.info(f"Found matching SIMSI summaries folder: {simsi_evidence_file}")

    if not os.path.isfile(simsi_evidence_file):
        simsi_evidence_file = extract_simsi_evidence_file_from_archive(simsi_evidence_file)
    
    return simsi_evidence_file


def extract_simsi_evidence_file_from_archive(simsi_evidence_file: Path) -> io.BytesIO:
    zip_file = Path(str(Path(simsi_evidence_file).parent.parent) + '.zip')
    if not zip_file.is_file():
        raise ValueError("Could not find SIMSI summaries folder nor archive.")
    
    logger.info(f'Extracting SIMSI evidence file from archive {zip_file}')
    archive = zipfile.ZipFile(zip_file, 'r')

    file_path_in_archive = zip_file.with_suffix('').name + '/p10/p10_evidence.txt'
    if not file_path_in_archive in archive.namelist():
        raise ValueError("Could not find SIMSI evidence file in archive.")
    return io.BytesIO(archive.read(file_path_in_archive))


def store_results_for_reuse(simsi_output_folder: Path, meta_input_file: Path, results_folder_name: str):
    """
    Rename the maracluster_output and summaries directories so we don't overwrite results
    """
    simsi_summaries_folder = simsi_output_folder / Path('summaries')
    simsi_summaries_folder_new = simsi_output_folder / Path(f'summaries_{results_folder_name}')
    simsi_summaries_folder.rename(simsi_summaries_folder_new)

    simsi_maracluster_folder = simsi_output_folder / Path('maracluster_output')
    simsi_maracluster_folder_new = simsi_output_folder / Path(f'maracluster_output_{results_folder_name}')
    simsi_maracluster_folder.rename(simsi_maracluster_folder_new)

    # keep a copy in the summaries folder so we can check in future runs if we can re-use the results
    copy_with_subprocess(str(meta_input_file), str(simsi_summaries_folder_new / meta_input_file.name))


def run_simsi_single(meta_input_file: Path,
                     simsi_output_folder: Path,
                     tmt_ms_level: str,
                     stringencies: int,
                     tmt_requantify: bool,
                     maximum_pep: int,
                     num_threads: int):
    # TODO: test that this works
    simsi_configs = ['--meta_input_file', str(meta_input_file),
                     '--output_folder', str(simsi_output_folder),
                     '--tmt_ms_level', str(tmt_ms_level),
                     '--tmt_requantify',
                     '--stringencies', str(stringencies),
                     '--maximum_pep', str(maximum_pep),
                     '--num_threads', str(num_threads)]
    if tmt_requantify:
        simsi_configs.append('--tmt_requantify')
    simsi.main(simsi_configs)


def create_meta_input_file(summary_files: List[str], data_type: str, simsi_folder: str, meta_input_file: Path):
    mq_txt_folders = [Path(f).parent for f in summary_files]
    experiments = []
    for summary in summary_files:
        summary_df = pd.read_csv(summary, sep='\t')
        experiments.append(summary_df['Experiment'][0])

    raw_folders = [get_raw_file_folder(simsi_folder, data_type, experiment) for experiment in experiments]

    correction_file_mapping_file = Path(simsi_folder) / 'correction_factor_files.tsv'

    # TODO: find better way to specify these paths
    mq_queue_file = Path(simsi_folder).parent / 'Queue' / 'queue.csv'
    correction_file_folder = Path(simsi_folder).parent / 'Queue' / 'Maxquant' / 'TMT_correction_factors'

    tmt_correction_files = get_correction_factor_files(experiments, correction_file_mapping_file, mq_queue_file,
                                                       correction_file_folder)

    meta_input_file_df = pd.DataFrame({'mq_txt_folder': mq_txt_folders,
                                       'raw_folder': raw_folders,
                                       'tmt_correction_file': tmt_correction_files})
    # meta_input_file_df = pd.DataFrame({'mq_txt_folder': mq_txt_folders,
    #                                    'raw_folder': raw_folders})
    meta_input_file_df.to_csv(meta_input_file, sep='\t', index=False)


def get_summary_files(raw_data_location, data_type, sample_annotation_file):
    logger.info(f"Getting paths to summary files {data_type}")
    try:
        sample_annotation_df = pd.read_csv(sample_annotation_file)
    except PermissionError:
        raise PermissionError("Cannot open patient annotation file, check if you have it open in Excel.")
    batches = sample_annotation.get_unique_batches(sample_annotation_df)

    summary_files = prep.get_data_location(raw_data_location, data_type, file_type='summary.txt')

    summary_files = prep.filter_evidence_files(summary_files, data_type.upper(), batches)

    return summary_files


def copy_files(summary_files, data_type, simsi_folder):
    """
    Copies raw files (skips already copied files) for all batches listed in sample_annotion based on the summary.txt of the MaxQuant search.
    """
    for summary in summary_files:
        logger.info(f"Copying raw files for {summary}")
        summary_df = pd.read_csv(summary, sep='\t')

        raw_files = summary_df['Raw file']
        raw_files = raw_files.loc[raw_files != 'Total']

        experiment = summary_df['Experiment'][0]

        batch_simsi_folder = get_raw_file_folder(simsi_folder, data_type, experiment)
        if not os.path.exists(batch_simsi_folder):
            os.makedirs(batch_simsi_folder)

        machines = ['lumos_1', 'eclipse_1', 'lumos_3', 'eclipse_2', 'lumos_2']
        if data_type == 'pp':
            machines = ['lumos_2', 'eclipse_2']

        for f in raw_files:
            # if raw file already exists as mzML don't copy
            if not os.path.exists(os.path.join(simsi_folder, 'simsi_output', data_type.upper(), 'mzML', f'{f}.mzML')):
                for machine in machines:
                    success = copy_with_subprocess(os.path.join('/media', 'raw_files', machine, 'raw', f'{f}.raw'), os.path.join(
                        batch_simsi_folder, f'{f}.raw'))
                    if success:
                        break
                if not success:
                    raise RuntimeError("Copying raw files for SIMSI failed")


def get_raw_file_folder(simsi_folder, data_type, experiment):
    return os.path.join(simsi_folder, 'raw_files', data_type.upper(), experiment)


def copy_with_subprocess(source, dest):
    if os.path.exists(dest):
        # logger.info(f"Skipping dest file {dest}, already copied")
        return True

    if not os.path.exists(source):
        logger.info(f"Could not find source file {source}")
        return False

    logger.info(" ".join(['cp', source, dest]))
    proc = subprocess.Popen(['cp', source, dest], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc.wait() == 0  # wait() returns the process' exit code, 0 means success


"""
python3 -m bin.simsi -c config_patients.json
"""
if __name__ == "__main__":
    main(sys.argv[1:])

