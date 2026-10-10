#!/usr/bin/env python3
"""Hash the source-only assessment and verify prior specialist seals unchanged."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "sleep_unified_research_v4"
OUT = V4 / "reviews/genomicsem_realistic_ld_feasibility_receipt_v4.json"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    assert not OUT.exists(), "preserve previous seal"
    prior = {
        "genomicsem_specialist_review_receipt_v4.json": "2fda9c9d2df98dacb0f2ecd43ab2465448fb377e617e2014ea627ff05e455f7d",
        "genomicsem_native_calibration_preparation_receipt_v4.json": "bc10d90c9f3211f1e8f15a3b3e31d062e9415604cf81a51a86077f4d20f25b32",
        "genomicsem_native_calibration_preparation_receipt_v4_2.json": "515bf63e6ebd750ab88ac77850c7d8fadb7c7deb5dd7a33b0fa50f3f7b54b7be",
        "genomicsem_native_calibration_preparation_receipt_v4_3.json": "83fa21c52721365bbeddcdbc343f62fc82435d4c644b49d8abd7f6eafd14cc00",
    }
    for name, expected in prior.items():
        seal = V4 / "reviews" / name
        assert sha(seal) == expected
        j = json.loads(seal.read_text())
        for r in j.get("records", j.get("artifacts", [])):
            p = Path(r["path"])
            if not p.is_absolute():
                p = ROOT / p
            assert p.stat().st_size == r["bytes"] and sha(p) == r["sha256"], p
    native_plan = V4 / "manifests/native_canonical_calibration_plan_v4_3.json"
    assert sha(native_plan) == "cfe5ef24f3455cc0508785e69b291d023f524b2e4872754598b864365f884c21"
    independent = V4 / "reviews/independent_canonical_executor_correction_review_v4_3.md"
    assert sha(independent) == "d165419725f70825fdc99aead77042e91720c052f31c3154f2f08d229f72528f"
    files = list((V4 / "statistical_validation/canonical_ld_feasibility_sources_v4").iterdir())
    files += [V4 / p for p in [
        "reviews/genomicsem_realistic_ld_calibration_feasibility_v4.md",
        "manifests/canonical_realistic_ld_source_admission_plan_v4.json",
        "statistical_validation/canonical_realistic_ld_feasibility_probe_v4.py",
        "statistical_validation/canonical_realistic_ld_feasibility_probe_receipt_v4.json",
        "statistical_validation/canonical_signed_genotype_generator_fixture_v4.py",
        "statistical_validation/canonical_signed_genotype_generator_fixture_receipt_v4.json",
        "statistical_validation/canonical_realistic_ld_source_admission_plan_v4.py",
        "statistical_validation/canonical_realistic_ld_feasibility_gates_v4.tsv",
        "statistical_validation/canonical_realistic_ld_feasibility_seal_v4.py",
    ]]
    rows = [{"path": str(p), "bytes": p.stat().st_size, "sha256": sha(p)} for p in sorted(files)]
    j = {"schema": "genomicsem-realistic-ld-feasibility-seal-v4.1", "sealed_utc": datetime.now(timezone.utc).isoformat(),
         "artifacts": rows, "artifact_count": len(rows), "total_artifact_bytes": sum(r["bytes"] for r in rows),
         "previous_specialist_seals_sha256": prior, "all_previous_artifacts_verified_unchanged": True,
         "native_control_plan_sha256": sha(native_plan), "independent_native_operational_review_sha256": sha(independent),
         "source_candidates_characterized": True, "restricted_signed_factor_qualification_proposed": True,
         "restricted_signed_factor_qualification_admitted": False, "full_realistic_LD_calibration_admitted": False,
         "invented_genotype_algebra_fixture_checks": 41, "source_genotype_decodes": 0,
         "new_source_acquisitions": 0, "native_estimator_or_empirical_calibration_jobs": 0,
         "new_41_covariance_outcomes": 0, "calibrated_heterogeneity_P_values": 0,
         "manuscript_or_historical_result_changes": False}
    OUT.write_text(json.dumps(j, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"seal": str(OUT), "sha256": sha(OUT), "artifacts": len(rows), "bytes": j["total_artifact_bytes"], "admitted": False}))


if __name__ == "__main__":
    main()
