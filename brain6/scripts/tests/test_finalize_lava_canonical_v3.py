from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/finalize_lava_canonical_v3.py"
spec = importlib.util.spec_from_file_location("finalize_lava_canonical_v3", SCRIPT)
finalizer = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(finalizer)


def base_inputs(qc_pass: bool):
    audit = {
        "analysis_id": "test-analysis",
        "run_id": "test-run",
        "state": "COMPLETE",
        "planned_loci": 2,
        "verified_loci": 2,
        "planned_cells": 4,
        "verified_cells": 4,
        "invalid_receipt_count": 0,
        "family_qc_pass": qc_pass,
        "currently_uncompleted_or_untested_cells": 1 if not qc_pass else 0,
        "maximum_allowed_untested_cells": 0,
        "statuses": {"TESTED": 3 if not qc_pass else 4, "NOT_RUN": 1 if not qc_pass else 0, "FAILED": 0},
        "family_lock_sha256": "family-hash",
        "execution_lock_sha256": "execution-hash",
        "canonical_script_sha256": "script-hash",
    }
    family = {
        "canonical_univariate": {"n_tests": 4, "loci": 2},
        "univariate_gate": {
            "familywise_alpha": 0.05,
            "n_tests": 4,
            "p_threshold_strictly_less_than": 0.0125,
        },
    }
    identity = {"run_id": "test-run", "run_identity_sha256": "identity-hash",
                "canonical_manifest_sha256": "manifest-hash"}
    aggregate = {"path": "/tmp/canonical.tsv", "sha256": "aggregate-hash", "rows": 4}
    return audit, family, identity, aggregate


def test_complete_qc_pass_authorizes_pairwise_but_does_not_run_it():
    audit, family, identity, aggregate = base_inputs(qc_pass=True)
    decision = finalizer.build_decision(audit, family, identity, aggregate, aggregate_rows=4)
    assert decision["overall_status"] == "PASS"
    assert decision["promotion_permitted"] is True
    assert decision["pairwise_stage_authorized"] is True
    assert decision["pairwise_stage_executed"] is False
    assert decision["multiple_testing"]["strict_p_threshold"] == 0.0125


def test_complete_qc_failure_is_finalized_without_promoting_results():
    audit, family, identity, aggregate = base_inputs(qc_pass=False)
    decision = finalizer.build_decision(audit, family, identity, aggregate, aggregate_rows=4)
    assert decision["overall_status"] == "FAILED_QC_NOT_PROMOTED"
    assert decision["promotion_permitted"] is False
    assert decision["pairwise_stage_authorized"] is False
    assert decision["multiple_testing"]["inference_promoted"] is False


def test_partial_family_cannot_be_finalized():
    audit, family, identity, aggregate = base_inputs(qc_pass=False)
    audit["state"] = "PARTIAL"
    with pytest.raises(ValueError, match="incomplete canonical family"):
        finalizer.build_decision(audit, family, identity, aggregate, aggregate_rows=4)


def test_runner_identity_hash_verifies_distinct_audit_and_identity_document_schemas():
    audit, _family, _identity, _aggregate = base_inputs(qc_pass=False)
    audit.update({
        "worker_count": 4,
        "input_provenance_sha256": "input-hash",
        "reference_provenance_sha256": "reference-hash",
    })
    identity_doc = {
        "analysis_id": audit["analysis_id"],
        "family_lock_sha256": audit["family_lock_sha256"],
        "execution_lock_sha256": audit["execution_lock_sha256"],
        "script_sha256": audit["canonical_script_sha256"],
        "runtime_validator_sha256": "validator-hash",
        "runtime_validation_output_sha256": "runtime-output-hash",
        "runtime": {"R": "4.3.3", "LAVA": "0.1.5"},
        "canonical_manifest_sha256": "manifest-hash",
    }
    runtime_identity = {
        "analysis_id": audit["analysis_id"],
        "family_lock_sha256": audit["family_lock_sha256"],
        "execution_lock_sha256": audit["execution_lock_sha256"],
        "input_provenance_sha256": audit["input_provenance_sha256"],
        "reference_provenance_sha256": audit["reference_provenance_sha256"],
        "batch_worker_sha256": identity_doc["script_sha256"],
        "runtime_validator_sha256": identity_doc["runtime_validator_sha256"],
        "runtime_validation_output_sha256": identity_doc["runtime_validation_output_sha256"],
        "R": "4.3.3",
        "LAVA": "0.1.5",
        "coordinator_sha256": "coordinator-hash",
    }
    run_id = hashlib.sha256(finalizer.canonical.canonical_json(runtime_identity)).hexdigest()
    audit["run_id"] = run_id
    identity_doc.update({"run_id": run_id, "run_identity_sha256": run_id})
    latest = {
        "analysis_id": audit["analysis_id"],
        "run_id": run_id,
        "run_identity": runtime_identity,
        "manifest_sha256": identity_doc["canonical_manifest_sha256"],
        "workers": 4,
    }

    # The runner records full identity inputs; its separate identity file is a
    # projection, so structural equality is neither expected nor appropriate.
    finalizer.verify_runner_identity(latest, identity_doc, audit)

    latest["run_identity"] = {**runtime_identity, "coordinator_sha256": "tampered"}
    with pytest.raises(ValueError, match="different run identity"):
        finalizer.verify_runner_identity(latest, identity_doc, audit)


def test_aggregate_rows_match_receipt_bound_cells_and_schema(tmp_path):
    run_dir = tmp_path / "run"
    (run_dir / "cells").mkdir(parents=True)
    cell_fields = ["phen", "locus_id", "status", "n_snps", "n_components", "h2.obs", "h2.latent", "p", "reason"]
    cells = [
        {"phen": trait, "locus_id": "1", "status": "TESTED", "n_snps": "50", "n_components": "2",
         "h2.obs": "0.1", "h2.latent": "0.1", "p": "0.01", "reason": ""}
        for trait in ("a", "b")
    ]
    with (run_dir / "cells/locus_1.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cells)
    aggregate = tmp_path / "aggregate.tsv"
    rows = [
        {"phen": row["phen"], "locus_id": "1", "chromosome": "1",
         **{field: row[field] for field in finalizer.CELL_FIELDS}}
        for row in cells
    ]
    with aggregate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=finalizer.AGGREGATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    assert finalizer.verify_aggregate(aggregate, run_dir, [{"LOC": "1", "CHR": "1"}], ("a", "b")) == 2

    rows[1]["p"] = "0.02"
    with aggregate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=finalizer.AGGREGATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="differs from receipt-bound cells"):
        finalizer.verify_aggregate(aggregate, run_dir, [{"LOC": "1", "CHR": "1"}], ("a", "b"))


def test_decision_write_is_immutable_and_idempotent(tmp_path):
    path = tmp_path / "decision.json"
    payload = finalizer.json_bytes({"status": "FAILED_QC_NOT_PROMOTED"})
    first = finalizer.write_immutable(path, payload)
    assert finalizer.write_immutable(path, payload) == first
    with pytest.raises(ValueError, match="Refusing to overwrite"):
        finalizer.write_immutable(path, finalizer.json_bytes({"status": "PASS"}))
    assert json.loads(path.read_text()) == {"status": "FAILED_QC_NOT_PROMOTED"}


def test_incomplete_live_run_creates_no_decision(tmp_path):
    with patch.object(finalizer.partial, "audit", return_value={"state": "PARTIAL"}), pytest.raises(
        ValueError, match="before all frozen cells are audited"
    ):
        finalizer.finalize(tmp_path, tmp_path, tmp_path / "reference.json")
    assert not (tmp_path / "canonical_family_decision.json").exists()
