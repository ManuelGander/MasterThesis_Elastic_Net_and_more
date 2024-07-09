import pandas as pd
import numpy as np

from bin import clinical_tools, config


def test_phospho_annot():
    results_folder = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.12.11_FH_minimal"
    config_file = f"{results_folder}/configs.json"
    preprocessed_pp_file = f"{results_folder}/preprocessed_pp.csv"
    reference_result_file = (
        f"{results_folder}/integration_tests/preprocessed_pp_phospho_annot.csv"
    )

    index_col = "Modified sequence"
    keep_default_na = False
    configs = config.load(config_file)
    preprocessed_pp = pd.read_csv(
        preprocessed_pp_file, index_col=index_col, keep_default_na=keep_default_na
    )

    patient_cols = preprocessed_pp.filter(regex=configs["patient_regex"]).columns
    preprocessed_pp[patient_cols] = (
        preprocessed_pp[patient_cols].replace("", np.nan).astype("float")
    )
    preprocessed_pp = clinical_tools.phospho_annot(
        preprocessed_pp,
        pspFastaFile=configs["clinic_proc"]["pspFastaFile"],
        pspKinaseSubstrateFile=configs["clinic_proc"]["pspKinaseSubstrateFile"],
        pspAnnotationFile=configs["clinic_proc"]["pspAnnotationFile"],
        pspRegulatoryFile=configs["clinic_proc"]["pspRegulatoryFile"],
    )

    # uncomment this to create new reference file
    preprocessed_pp.to_csv(reference_result_file, sep="\t", index=False)

    preprocessed_pp_reference = pd.read_csv(
        reference_result_file,
        sep="\t",
        dtype={"PSP_LT_LIT": "object", "PSP_MS_LIT": "object", "PSP_MS_CST": "object"},
    )
    non_patient_cols = [
        c for c in preprocessed_pp_reference.columns if c not in patient_cols
    ]
    preprocessed_pp_reference[non_patient_cols] = preprocessed_pp_reference[
        non_patient_cols
    ].fillna("")
    
    pd.testing.assert_frame_equal(
        preprocessed_pp,
        preprocessed_pp_reference,
        check_like=True,
        check_dtype=False,
        check_exact=False,
    )


if __name__ == "__main__":
    test_phospho_annot()
