#!/usr/bin/env python3
"""Check protected legacy B method and row eligibility against frozen v3 rules."""
from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/track_b_admission_v1"
ARCHIVE = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
REL = Path("results/track_b/pleiotropy")
RUN = "f8a3de8642698a84cddaf59aa2d94e518a5994139fe8fc2bb5ff530ccdaf5066"
BUNDLE = ARCHIVE / REL / "continuations/post_lava_terminal_v2/sequential_execution_records" / RUN / "B/reproducibility_and_cleanup.plan.json"
LEDGER = ARCHIVE / REL / "results/placo/B.full.tsv.gz"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    criteria_path = OUT / "criteria.lock.json"
    criteria = json.loads(criteria_path.read_text())
    expected = criteria["expected"]
    preflight_path = OUT / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    if preflight["status"] != "PASS" or preflight["criteria_sha256"] != sha(criteria_path):
        raise ValueError("locked preflight is not passing")
    if sha(BUNDLE) != expected["cleanup_bundle_sha256"] or sha(LEDGER) != expected["ledger_sha256"]:
        raise ValueError("historical source hashes changed")
    bundle = json.loads(BUNDLE.read_text())
    info = bundle["embedded_reproducibility_artifacts"]["materialization.provenance.json"]
    raw = base64.b64decode(info["base64"], validate=True)
    if hashlib.sha256(raw).hexdigest() != info["sha256"] or len(raw) != info["bytes"]:
        raise ValueError("materialization receipt content changed")
    mat = json.loads(raw)
    sidecar = json.loads((ARCHIVE / REL / "results/placo/B.provenance.json").read_text())
    family = json.loads((ROOT / "brain6/config/placo_family_v3/family_lock.json").read_text())
    policy = json.loads((ROOT / "config/track_b_pleiotropy_policy.json").read_text())
    counts = mat["alignment_counts"]
    sources = {x["trait_id"]: x for x in mat["source_files"]}
    checks = {
        "exact_run_pair_and_source_receipt": mat["run_fingerprint"] == RUN and mat["pair_id"] == "B"
            and mat["analysis_id"] == "track-b-v1.0-pleiotropy"
            and sidecar["materialization_provenance_sha256"] == info["sha256"]
            and mat["aligned_output"]["sha256"] == sidecar["input_sha256"],
        "locked_gwas_and_eligibility": mat["trait1"] == "insomnia" and mat["trait2"] == "adhd"
            and set(sources) == {"insomnia", "adhd"}
            and all(sources[t]["sha256"] == expected[f"{t}_dense_sha256"] for t in sources)
            and mat["aligned_output"]["rows"] == counts["eligible_written"] == expected["rows"]
            and counts["allele_matches"] == counts["eligible_written"] + counts["z_squared_exclusions"]
            and counts["invalid_z_or_p_dropped"] == counts["duplicate_ids_dropped"] == 0
            and len(mat["per_autosome"]) == 22
            and all(x["eligible_written"] > 0 for x in mat["per_autosome"]),
        "official_correlated_placo_nuisance": mat["placo_source_sha256"] == family["placo_source_sha256"] == expected["placo_source_sha256"]
            and sidecar["nuisance_estimation"]["method"] == "PINNED_PLACO_PLUS_var.placo_AND_cor.pearson"
            and family["method"] == "PLACO_PLUS"
            and family["p_threshold"] == policy["placo_plus"]["marginal_p_threshold_for_nuisance_estimation"] == criteria["frozen_statistics"]["marginal_p_threshold"],
        "historical_tolerance_retained_as_legacy": policy["placo_plus"]["absolute_tolerance"] == criteria["frozen_statistics"]["historical_absolute_tolerance"]
            and family["absolute_tolerance"] == criteria["frozen_statistics"]["v3_absolute_tolerance"]
            and criteria["frozen_statistics"]["historical_abs_tol_difference_policy"]
                == "DOCUMENT_AS_PROTECTED_LEGACY_METHOD_PARAMETER;DO_NOT_RECOMPUTE_OR_REPLACE_RESULTS",
        "five_track_rule_frozen": family["n_selected_tracks"] == criteria["frozen_statistics"]["n_selected_tracks"] == 5
            and family["headline_threshold"] == criteria["frozen_statistics"]["headline_p_strictly_less_than"]
            and family["extreme_z2"] == policy["placo_plus"]["z_squared_maximum"] == criteria["frozen_statistics"]["z_squared_maximum"]
            and family["protected_legacy_pair"] == "insomnia__adhd",
    }
    # Audit the preserved row family independently of the original materializer's summary.
    max_z2 = 0.0
    rows = 0
    with gzip.open(LEDGER, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            rows += 1
            z2 = max(float(row["Z1"]) ** 2, float(row["Z2"]) ** 2)
            if not math.isfinite(z2):
                raise ValueError(f"nonfinite Z in archived result at row {rows}")
            max_z2 = max(max_z2, z2)
    checks["every_retained_row_meets_extreme_z_rule"] = rows == expected["rows"] and max_z2 <= family["extreme_z2"]
    result = {"status":"PASS" if all(checks.values()) else "FAIL",
              "criteria_sha256":sha(criteria_path),"preflight_sha256":sha(preflight_path),
              "validator_sha256":sha(Path(__file__)),"historical_ledger_sha256":sha(LEDGER),
              "materialization_provenance_sha256":info["sha256"],"checks":checks,
              "failed_checks":[k for k,v in checks.items() if not v],
              "retained_rows":rows,"largest_retained_z_squared":max_z2,
              "historical_abs_tol":policy["placo_plus"]["absolute_tolerance"],
              "v3_abs_tol":family["absolute_tolerance"],
              "interpretation":"The protected legacy B result uses the same official correlated PLACO+ implementation and nuisance rules; its historical absolute tolerance is preserved as an explicit legacy parameter, and its stricter archived row family remains separate from the four v3 pair outputs."}
    receipt = OUT / "scientific_compatibility.json"
    payload = json.dumps(result,indent=2,sort_keys=True)+"\n"
    if receipt.exists() and receipt.read_text() != payload:
        raise FileExistsError("refusing to overwrite scientific compatibility receipt")
    if not receipt.exists():
        receipt.write_text(payload)
    print(json.dumps({"status":result["status"],"failed_checks":result["failed_checks"],"largest_retained_z_squared":max_z2}))
    if result["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
