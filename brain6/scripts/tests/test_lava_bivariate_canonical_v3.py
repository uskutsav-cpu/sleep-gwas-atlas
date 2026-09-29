from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/run_lava_bivariate_canonical_v3.py"
spec = importlib.util.spec_from_file_location("lava_bivariate_canonical_v3", SCRIPT)
runner = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(runner)


def make_family(*, failed: int = 0) -> list[dict[str, str]]:
    rows = []
    for index in range(12_475):
        status = "FAILED" if index < failed else "UNIVARIATE_UNDERPOWERED"
        p, rg = "NA", "NA"
        if index == failed:
            status, p, rg = "TESTED", "1e-6", "0.2"
        rows.append({"pair_id": runner.PAIRS[index // 2_495],
            "locus_id": str(index % 2_495 + 1), "status": status,
            "p": p, "local_rg": rg, "reason": ""})
    return rows


def test_family_q_uses_full_fixed_denominator_and_keeps_nontested_q_na():
    summary = runner.bivariate_summary(make_family(), 0.01)
    assert summary["planned_slots"] == summary["bh_family_size"] == 12_475
    assert summary["status"] == "PASS"
    tested = [row for row in summary["rows"] if row["status"] == "TESTED"]
    untested = [row for row in summary["rows"] if row["status"] != "TESTED"]
    assert tested[0]["family_fdr"] == pytest.approx(0.012475)
    assert all(row["family_fdr"] == "NA" for row in untested)


def test_execution_failure_ceiling_is_strict_and_preserves_all_slots():
    accepted = runner.bivariate_summary(make_family(failed=124), 0.01)
    rejected = runner.bivariate_summary(make_family(failed=125), 0.01)
    assert accepted["planned_slots"] == rejected["planned_slots"] == 12_475
    assert accepted["status"] == "INSUFFICIENT_EVIDENCE"
    assert accepted["execution_failure_fraction"] <= 0.01
    assert rejected["status"] == "FAILED_QC_NOT_CONSUMED"
    assert rejected["execution_failure_fraction"] > 0.01


def test_no_overlap_is_reported_but_does_not_consume_execution_failure_allowance():
    rows = make_family()
    for row in rows[1:126]:
        row["status"] = "NO_OVERLAP"
    summary = runner.bivariate_summary(rows, 0.01)
    assert summary["counts"]["NO_OVERLAP"] == 125
    assert summary["execution_failures_or_missing"] == 0
    assert summary["execution_failure_fraction"] == 0
    assert summary["status"] == "INSUFFICIENT_EVIDENCE"
    tested = next(row for row in summary["rows"] if row["status"] == "TESTED")
    assert tested["family_fdr"] == pytest.approx(0.012475)


def test_bivariate_lock_binds_canonical_family_worker_and_four_workers():
    lock, lock_sha = runner.validate_lock(runner.LOCK)
    canonical = json.loads((ROOT / "brain6/config/lava_family_canonical_v3.json").read_text())
    assert lock_sha == runner.sha256(runner.LOCK)
    assert lock["canonical_source"]["family_lock_sha256"] == runner.sha256(
        ROOT / "brain6/config/lava_family_canonical_v3.json")
    assert lock["canonical_source"]["maximum_untested_fraction"] == 0.05
    assert lock["pair_family"]["n_slots"] == 12_475
    assert lock["execution"]["worker_count"] == 4
    assert lock["pair_family"]["pair_ids"] == canonical["pairs"]


def test_canonical_rows_for_locus_requires_exact_seven_trait_coverage(tmp_path: Path):
    path = tmp_path / "cells.tsv"
    rows = [{"phen": trait, "locus_id": "9", "status": "TESTED", "p": "0.1"}
            for trait in runner.TRAITS]
    path.write_text("phen\tlocus_id\tstatus\tp\n" + "".join(
        f"{row['phen']}\t{row['locus_id']}\t{row['status']}\t{row['p']}\n" for row in rows))
    assert len(runner.canonical_rows_for_locus(path, "9")) == 7
    with pytest.raises(ValueError, match="exactly seven"):
        runner.canonical_rows_for_locus(path, "10")


def test_canonical_aggregate_is_compared_row_for_row_to_receipt_cells(tmp_path: Path):
    canonical_root = tmp_path / "canonical"
    cells_dir = canonical_root / "cells"
    cells_dir.mkdir(parents=True)
    cell_fields = ["phen", "locus_id", "status", "n_snps", "n_components",
                   "h2.obs", "h2.latent", "p", "reason"]
    cells = [{"phen": trait, "locus_id": "9", "status": "TESTED", "n_snps": "20",
              "n_components": "3", "h2.obs": "0.1", "h2.latent": "0.2", "p": "0.01", "reason": ""}
             for trait in runner.TRAITS]
    with (cells_dir / "locus_9.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(cells)
    aggregate = tmp_path / "canonical_family_results.tsv"
    fields = list(runner.CANONICAL_AGGREGATE_FIELDS)
    expected = [{**row, "chromosome": "1"} for row in cells]
    with aggregate.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(expected)
    loci = [{"LOC": "9", "CHR": "1", "START": "10", "STOP": "20"}]
    runner.verify_canonical_aggregate(aggregate, canonical_root, loci)

    expected[0]["p"] = "0.02"
    with aggregate.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(expected)
    with pytest.raises(ValueError, match="differs from verified cells"):
        runner.verify_canonical_aggregate(aggregate, canonical_root, loci)


def test_pair_coverage_checks_exact_threshold_and_inventory():
    rows = [
        {"pair_id": runner.PAIRS[0], "locus_id": "9", "shared_reference_variants": "2", "status": "READY"},
        {"pair_id": runner.PAIRS[1], "locus_id": "9", "shared_reference_variants": "1", "status": "NO_OVERLAP_LT_MIN_K"},
    ]
    expected = {(runner.PAIRS[0], "9"), (runner.PAIRS[1], "9")}
    assert runner.validate_pair_coverage(rows, expected, 2) == {
        (runner.PAIRS[0], "9"): 2, (runner.PAIRS[1], "9"): 1}
    with pytest.raises(ValueError, match="status disagrees"):
        runner.validate_pair_coverage([dict(rows[1], status="READY")],
            {(runner.PAIRS[1], "9")}, 2)
    with pytest.raises(ValueError, match="Negative"):
        runner.validate_pair_coverage([dict(rows[1], shared_reference_variants="-1")],
            {(runner.PAIRS[1], "9")}, 2)
    with pytest.raises(ValueError, match="exactly cover"):
        runner.validate_pair_coverage(rows[:1], expected, 2)


def test_pair_gate_outcomes_obey_strict_gate_before_r_startup():
    lock, _ = runner.validate_lock(runner.LOCK)
    threshold = lock["eligibility"]["strict_p_threshold"]
    rows = [{"phen": trait, "locus_id": "9", "status": "TESTED", "p": "1e-7",
             "h2.obs": "0.1", "h2.latent": "0.2", "reason": ""} for trait in runner.TRAITS]
    by_trait = {row["phen"]: row for row in rows}
    by_trait["longsleep"]["p"] = str(threshold)
    pairs = [{"pair_id": pair, "locus_id": "9", "phenotypes": pair.split("__"),
              "shared_reference_variants": 2} for pair in runner.PAIRS]
    gated, eligible = runner.pair_gate_outcomes(pairs, rows, 2, threshold)
    assert eligible == ["insomnia__adhd", "insomnia__mdd"]
    assert {row["pair_id"] for row in gated} == {"longsleep__scz", "longsleep__bipolar", "longsleep__parkinson"}
    assert all(row["status"] == "UNIVARIATE_UNDERPOWERED" for row in gated)
    by_trait["insomnia"]["p"] = str(threshold)
    gated, eligible = runner.pair_gate_outcomes(pairs, rows, 2, threshold)
    assert eligible == []
    assert len(gated) == 5
    pairs[0]["shared_reference_variants"] = 1
    gated, eligible = runner.pair_gate_outcomes(pairs, rows, 2, threshold)
    assert "insomnia__adhd" not in eligible
    assert next(row for row in gated if row["pair_id"] == "insomnia__adhd")["status"] == "NO_OVERLAP"
    malformed_pairs = [dict(pair) for pair in pairs]
    malformed_pairs[0]["phenotypes"] = ["adhd", "insomnia"]
    with pytest.raises(ValueError, match="identifier does not match"):
        runner.pair_gate_outcomes(malformed_pairs, rows, 2, threshold)
    malformed_rows = [dict(row) for row in rows]
    malformed_rows[0]["status"] = "UNRECOGNIZED"
    with pytest.raises(ValueError, match="unknown canonical status"):
        runner.pair_gate_outcomes(pairs, malformed_rows, 2, threshold)


def test_gate_only_locus_receipt_is_complete_and_resumable_without_r(tmp_path: Path):
    canonical_rows = [{"phen": trait, "locus_id": "9", "status": "NOT_RUN",
        "h2.obs": "NA", "h2.latent": "NA", "p": "NA", "reason": "no canonical result"}
        for trait in runner.TRAITS]
    pair_results = [{"pair_id": pair, "locus_id": "9", "status": "UNIVARIATE_UNDERPOWERED",
        "p": "NA", "local_rg": "NA", "reason": "ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_NOT_RUN"}
        for pair in runner.PAIRS]
    phenotype_pairs = {"insomnia__adhd": ("insomnia", "adhd"),
        "insomnia__mdd": ("insomnia", "mdd"), "longsleep__scz": ("longsleep", "scz"),
        "longsleep__bipolar": ("longsleep", "bipolar"),
        "longsleep__parkinson": ("longsleep", "parkinson")}
    pair_gate_inputs = [{"pair_id": pair, "locus_id": "9", "phenotypes": list(phenotype_pairs[pair]),
        "shared_reference_variants": 10} for pair in runner.PAIRS]
    target = tmp_path / "loci" / "9"
    outcome = runner.write_gate_only_unit(target, "9", "b" * 64, "c" * 64, "d" * 64,
        "e" * 64, "f" * 64, canonical_rows, pair_results, pair_gate_inputs,
        2, 2.86286859433152e-06, "a" * 64, 2)
    assert outcome["execution_mode"] == "PYTHON_GATE_ONLY"
    assert outcome["attempts"] == 0
    receipt = runner.verify_unit(target, "9", "b" * 64, "c" * 64, "c" * 64, "d" * 64, "e" * 64)
    assert receipt and receipt["execution_mode"] == "PYTHON_GATE_ONLY"
    config = json.loads((target / "worker_config.json").read_text())
    assert config["r_process_attempts"] == 0


def test_coordinator_does_not_start_r_when_every_pair_is_gated(tmp_path: Path, monkeypatch):
    locus_id = "9"
    canonical_root = tmp_path / "canonical"
    canonical_root.mkdir()
    cells_path = runner.canonical.result_path(canonical_root, locus_id)
    cell_fields = list(runner.CANONICAL_AGGREGATE_FIELDS)
    cells = [{"phen": trait, "locus_id": locus_id, "chromosome": "1", "status": "NOT_RUN",
        "n_snps": "NA", "n_components": "NA", "h2.obs": "NA", "h2.latent": "NA",
        "p": "NA", "reason": "MISSING_TEST"} for trait in runner.TRAITS]
    cells_path.parent.mkdir(parents=True)
    with cells_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cell_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(cells)
    canonical_receipt_path = runner.canonical.receipt_path(canonical_root, locus_id)
    runner.atomic_json(canonical_receipt_path, {"output": {"sha256": runner.sha256(cells_path)}})
    lock, lock_sha = runner.validate_lock(runner.LOCK)
    family = json.loads((ROOT / "brain6/config/lava_family_canonical_v3.json").read_text())
    pair_lock = {pair["pair_id"]: pair for pair in json.loads(
        (ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json").read_text())["pairs"]}
    source = {"run_id": "d" * 64,
        "coverage": {(pair, locus_id): 10 for pair in runner.PAIRS}}
    def unexpected_r(*_args, **_kwargs):
        raise AssertionError("R must not start when no pair passes the frozen gate")
    monkeypatch.setattr(runner.subprocess, "run", unexpected_r)
    outcome = runner.write_worker_unit(tmp_path / "bivariate" / locus_id,
        {"LOC": locus_id, "CHR": "1"}, family, lock, lock_sha, "b" * 64,
        canonical_root, source, tmp_path / "inputs", tmp_path / "reference",
        pair_lock, Path("Rscript"), "runtime")
    assert outcome["execution_mode"] == "PYTHON_GATE_ONLY"
    assert outcome["attempts"] == 0
    assert (tmp_path / "bivariate" / locus_id / "receipt.json").is_file()
    resumed = runner.write_worker_unit(tmp_path / "bivariate" / locus_id,
        {"LOC": locus_id, "CHR": "1"}, family, lock, lock_sha, "b" * 64,
        canonical_root, source, tmp_path / "inputs", tmp_path / "reference",
        pair_lock, Path("Rscript"), "runtime")
    assert resumed["resumed"] is True
    assert resumed["execution_mode"] == "PYTHON_GATE_ONLY"


def test_bivariate_stage_fails_closed_without_a_completed_canonical_audit(tmp_path: Path):
    lock, _ = runner.validate_lock(runner.LOCK)
    canonical_root = tmp_path / "canonical"
    canonical_root.mkdir()
    source = lock["canonical_source"]
    (canonical_root / "run_identity.json").write_text(json.dumps({
        "run_id": source["run_id"], "run_identity_sha256": source["run_id"],
        "family_lock_sha256": source["family_lock_sha256"],
        "execution_lock_sha256": source["execution_lock_sha256"],
    }))
    with pytest.raises(RuntimeError, match="no completed full-family audit"):
        runner.verify_canonical_source(lock, canonical_root,
            ROOT / "brain6/config/lava_family_canonical_v3.json",
            ROOT / "brain6/config/lava_execution_canonical_v3.json",
            tmp_path / "inputs", tmp_path / "reference")
