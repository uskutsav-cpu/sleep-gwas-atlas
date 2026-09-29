from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "build_grover_2022_sleep_pd_ldsc_context.py"
SPEC = importlib.util.spec_from_file_location("build_grover_2022_sleep_pd_ldsc_context", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_supplement_identity_is_hash_bound_without_claiming_contents_were_reviewed():
    row, provenance = MODULE.build()
    source = provenance["source"]
    supplement = source["linked_supplementary_file"]
    assert row["numeric_pair_estimate_reviewed"] == "False"
    assert source["linked_supplementary_workbook_reviewed"] is False
    assert supplement["filename"] == "Table_1.xlsx"
    assert supplement["mime_type"] == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    assert supplement["size_bytes"] == 722516
    assert supplement["binary_retrieved"] is False
    assert supplement["binary_sha256"] is None
    assert supplement["numeric_pair_estimate_reviewed"] is False
    assert provenance["source_file_hashes"][str(MODULE.ARTICLE_XML.relative_to(MODULE.ROOT))] == MODULE.sha256(MODULE.ARTICLE_XML)
