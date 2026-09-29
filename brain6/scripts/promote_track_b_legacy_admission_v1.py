#!/usr/bin/env python3
"""Promote exactly one validated historical B result into the protected Brain6 slot."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/track_b_admission_v1"
ARCHIVE = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
SOURCE = ARCHIVE / "results/track_b/pleiotropy/results/placo"
FINAL = ROOT / "brain6/results/placo/insomnia__adhd"
TEMP = FINAL.parent / ".insomnia__adhd.legacy_admission_v1.staging"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def require(predicate: bool, message: str) -> None:
    if not predicate:
        raise ValueError(message)


def main() -> None:
    criteria_path = OUT / "criteria.lock.json"
    criteria = read(criteria_path)
    expected = criteria["expected"]
    preflight_path = OUT / "preflight.json"
    adapter_path = OUT / "adapter_validation.json"
    scientific_path = OUT / "scientific_compatibility.json"
    preflight, adapter, scientific = map(read, (preflight_path, adapter_path, scientific_path))
    criteria_hash = sha(criteria_path)
    require(criteria["promotion_rule"] == "ALL_REQUIRED_CHECKS_TRUE;NO_CANONICAL_MUTATION_BEFORE_VALIDATOR_AND_ADAPTER_PASS", "promotion rule changed")
    require(preflight["status"] == adapter["status"] == scientific["status"] == "PASS", "required validation is not passing")
    require(not preflight["failed_checks"] and not scientific["failed_checks"] and all(x["pass"] for x in preflight["checks"].values()), "a required check failed")
    require(preflight["criteria_sha256"] == adapter["criteria_sha256"] == scientific["criteria_sha256"] == criteria_hash, "criteria identity drifted")
    require(preflight["validator_sha256"] == sha(ROOT / "brain6/scripts/validate_track_b_legacy_admission_v1.py"), "preflight validator changed")
    require(adapter["adapter_sha256"] == sha(ROOT / "brain6/scripts/adapter_track_b_legacy_v1.py"), "adapter changed")
    require(scientific["validator_sha256"] == sha(ROOT / "brain6/scripts/validate_track_b_scientific_compatibility_v1.py"), "scientific validator changed")
    require(adapter["preflight_sha256"] == scientific["preflight_sha256"] == sha(preflight_path), "preflight receipt changed")
    require(adapter["source_ledger_sha256"] == scientific["historical_ledger_sha256"] == expected["ledger_sha256"], "source ledger identity changed")
    require(adapter["rows"] == adapter["unique_snps"] == adapter["independent_bh_checked_rows"] == scientific["retained_rows"] == expected["rows"], "full-row validation incomplete")
    require(adapter["independent_bh_q_mismatches"] == adapter["five_track_adjustment_mismatches"] == 0, "multiplicity audit failed")
    require(adapter["source_fields"] == ["analysis_id","pair_id","SNP","CHR","BP","A1","A2","Z1","Z2","P1","P2","T_PLACO_PLUS","P_PLACO_PLUS","PLACO_BH_Q","within_pair_family_n","analysis_status","numerical_error"], "17-field source schema changed")
    require(adapter["adapted_fields"] == ["SNP","CHR","BP","Z1","Z2","P_PLACO","Q_WITHIN_PAIR","Q_BONFERRONI_ACROSS_PAIRS","status","headline"], "10-field adapted schema changed")
    require(adapter["pair_id"] == criteria["pair_id_brain6"] == "insomnia__adhd", "pair identity changed")
    view = Path(adapter["staged_view_path"])
    require(view.is_file() and sha(view) == adapter["staged_view_sha256"], "staged view hash changed")
    original = SOURCE / "B.full.tsv.gz"
    provenance = SOURCE / "B.provenance.json"
    require(sha(original) == expected["ledger_sha256"] and sha(provenance) == expected["provenance_sha256"], "historical bytes changed")
    family_path = ROOT / "brain6/config/placo_family_v3/family_lock.json"
    family = read(family_path)
    require(sha(family_path) == expected["family_lock_sha256"] and family["protected_legacy_pair"] == "insomnia__adhd", "family lock changed")
    four = read(ROOT / "brain6/results/placo/placo_v3_pair_qc_validation.json")
    require(four["family_lock_sha256"] == expected["family_lock_sha256"] and len(four["published_output_sha256"]) == 4, "four-pair QC receipt changed")
    for pair, digest in four["published_output_sha256"].items():
        require(sha(ROOT / "brain6/results/placo" / pair / "variants.tsv.gz") == digest, f"frozen v3 output changed: {pair}")
    require(not FINAL.exists() and not TEMP.exists(), "protected slot or staging directory already exists")
    TEMP.mkdir(parents=False, exist_ok=False)
    copies = ((original, "B.full.tsv.gz", expected["ledger_sha256"]),
              (provenance, "B.provenance.json", expected["provenance_sha256"]),
              (view, "variants.tsv.gz", adapter["staged_view_sha256"]))
    for source, name, digest in copies:
        destination = TEMP / name
        with source.open("rb") as incoming, destination.open("xb") as outgoing:
            shutil.copyfileobj(incoming, outgoing, length=8 << 20)
            outgoing.flush()
            os.fsync(outgoing.fileno())
        require(sha(destination) == digest, f"staged protected copy differs: {name}")
    receipt = {
        "analysis_id":criteria["analysis_id"], "status":"PROMOTED_AFTER_ALL_CHECKS_PASS",
        "admitted_at_utc":datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "historical_run_fingerprint":criteria["historical_run_fingerprint"],
        "pair_id":"insomnia__adhd","criteria_sha256":criteria_hash,
        "preflight_sha256":sha(preflight_path),"adapter_validation_sha256":sha(adapter_path),
        "scientific_compatibility_sha256":sha(scientific_path),
        "family_lock_sha256":expected["family_lock_sha256"],
        "files":{name:{"sha256":digest,"bytes":(TEMP / name).stat().st_size} for _,name,digest in copies},
        "historical_method_parameter_policy":criteria["frozen_statistics"]["historical_abs_tol_difference_policy"],
        "frozen_four_pair_output_sha256":four["published_output_sha256"],
        "interpretation":"Exact archived B ledger and provenance are retained beside the validated 10-field v3 view; the four prior v3 tracks were not rewritten."
    }
    (TEMP / "admission_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    require(not FINAL.exists(), "protected slot appeared during staging")
    TEMP.rename(FINAL)
    print(json.dumps({"status":receipt["status"],"pair_id":receipt["pair_id"],"protected_slot":str(FINAL),"rows":adapter["rows"]}))


if __name__ == "__main__":
    main()
