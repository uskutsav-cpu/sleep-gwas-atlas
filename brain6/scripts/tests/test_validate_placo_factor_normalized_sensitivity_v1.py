from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "validate_placo_factor_normalized_sensitivity_v1.py"
SPEC = importlib.util.spec_from_file_location("placo_factor_normalized_validator", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(module)


def edge(**overrides):
    row = {
        "CHR": "17", "SNP1": "rs1", "SNP2": "rs2", "lava_r_unscaled": "1.0001",
        "factor_r_normalized": "0.9999", "factor_r2_normalized": "0.99980001",
        "diagnostic_status": "NORMALIZED_FACTOR_SENSITIVITY_ONLY_NOT_GENOTYPE_VALIDATED",
    }
    return {**row, **overrides}


def locus(pair, chrom, start, stop, lead):
    return {
        "pair_id": pair, "CHR": str(chrom), "START": str(start), "STOP": str(stop),
        "lead_variant": lead, "lead_variants": lead,
    }


def test_range_exception_rows_require_out_of_range_raw_and_bounded_normalized_ld():
    assert module.validate_range_edges([edge()]) == 1
    with pytest.raises(ValueError, match="raw out-of-range"):
        module.validate_range_edges([edge(lava_r_unscaled="1.0")])
    with pytest.raises(ValueError, match=r"outside \[-1,1\]"):
        module.validate_range_edges([edge(factor_r_normalized="1.0000001")])


def test_range_exception_rows_require_explicit_diagnostic_only_status():
    with pytest.raises(ValueError, match="unexpected promotion"):
        module.validate_range_edges([edge(diagnostic_status="PROMOTED")])


def test_factor_normalized_sensitivity_matches_baseline_loci_and_leads():
    baseline = [locus("pair", 17, 100, 200, "rs1")]
    result = module.validate_loci_match_baseline(baseline, baseline)
    assert result["interval_sets_match"] == result["lead_sets_match"] == 1
    changed = [locus("pair", 17, 100, 200, "rs2")]
    result = module.validate_loci_match_baseline(baseline, changed)
    assert result["interval_sets_match"] == 1
    assert result["lead_sets_match"] == 0


def test_legacy_v1_output_integrity_is_preserved_and_fail_closed(tmp_path):
    out = tmp_path / "brain6/results/loci/placo_factor_normalized_sensitivity_v1"
    out.mkdir(parents=True)
    records = []
    for name in ("variants", "loci", "replaced_edges", "summary"):
        path = out / f"{name}.tsv"
        payload = f"{name}\n".encode()
        path.write_bytes(payload)
        records.append({"name": name, "path": str(path.relative_to(tmp_path)),
                        "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()})
    (out / "provenance.json").write_text(json.dumps({
        "analysis_id": "brain6_placo_factor_normalized_sensitivity_v1",
        "status": "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED",
        "outputs": records,
    }))
    assert set(module.validate_legacy_v1_integrity(tmp_path)) == {
        "variants", "loci", "replaced_edges", "summary"
    }
    (out / "loci.tsv").write_text("changed\n")
    with pytest.raises(ValueError, match="archived factor-normalized v1 output changed"):
        module.validate_legacy_v1_integrity(tmp_path)
