# WP3 Sample Pipeline

Automated sample analysis pipeline for data processing of TMT and label-free data. 

Input parameters for running the pipeline can be set in configs.json file.

The software is dockerized and so the easiest way to run it is by .....

## Requirements and installation

When not using docker
OS: Ubuntu >=18
Python: ">=3.9, <3.11"

## Usage

### How to run pipeline with make and poetry

This does not work with older pandas versions, so use the poetry environment!
```
poetry shell
export PYTHON_KEYRING_BACKEND=keyring.backends.null.Keyring
poetry install
```

To run the pipeline, adjust the file paths in `configs.json` and run:
```
make all
```

### How to run pipeline with docker

Login and pull image
```
docker login gitlab.lrz.de:5005
docker pull gitlab.lrz.de:5005/proteomics/topas/wp3_sample_pipeline:latest
```

To run the pipeline, adjust the input parameters and file paths in `config_patients.json` and run:
```
make docker_all
```

For minimal test run (3 batches) run:

```
CONFIG_FILE=config_patients_minimal_test.json make docker_all
```


### Comparing results from different runs

Run the script `compare_result_folders.py` with folders to compare as input arguments. Example:
```
python ./compare_result_folders.py --folder_1 <dir> --folder_2 <dir> 
```

This will create a folder with results in pdf reports in the `Retrospective_study/basket_score_comparisons` folder with correlation 
plots of the basket scores as well as detected outliers in the analysis of big changes of protein/p-site ranks.


### Prepare MQ results for SIMSI-transfer

Run the script `raw_files_for_simsi.py` with data folder of all batches as input as well as which cohort ('sarcoma', 'chordoma' or 'all') is 
wanted. 

Example:
```
python ./raw_files_for_simsi.py --data_folder <dir> --cohort <str> 
```

This will create files for full and phosphoproteome to use for rsync to copy raw files to folder of project
`{raw_data_location}/rsync_include_fp_{cohort}_{date}.txt`
