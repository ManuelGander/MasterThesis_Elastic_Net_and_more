

import bin.TOPAS_drug_scoring as drug_scoring


def test_read_regulated_curves_file():
    drug_scoring.read_regulated_curves_file('A204', [], '/tmp')