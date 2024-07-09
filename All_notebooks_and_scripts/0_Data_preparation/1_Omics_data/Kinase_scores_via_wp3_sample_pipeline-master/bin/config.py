import json


def load(config_file):
    with open(config_file, "r") as inp:
        configs = json.load(inp)

    if "data_types" not in configs:
        configs["data_types"] = ["fp", "pp"]

    if "run_lfq" not in configs["preprocessing"]:
        configs["preprocessing"]["run_lfq"] = False

    if "metadata_annotation" not in configs:
        configs["metadata_annotation"] = (
            "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/"
            "MasterInform_Metadata_MTB_Portal_221129.xlsx"
        )

    if "num_threads" not in configs["simsi"]:
        configs["simsi"]["num_threads"] = 4

    if "fasta_file" not in configs["preprocessing"]:
        configs["preprocessing"]["fasta_file"] = (
            "/media/kusterlab/internal_projects/active/TOPAS/Databases/"
            "uniprot_proteome_up000005640_03112020_cdkn2a_isoforms.fasta"
        )
    
    if "patient_regex" not in configs:
        configs["patient_regex"] = "^\\S+-.+-\\S+"

    if "portal" not in configs:
        configs["portal"] = {"update": 0}

    return configs
