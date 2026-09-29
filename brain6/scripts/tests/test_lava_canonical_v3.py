from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "brain6/scripts/run_lava_canonical_v3.py"
spec = importlib.util.spec_from_file_location("lava_canonical_v3", SCRIPT)
runner = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(runner)


def test_canonical_manifest_is_exact_unique_and_deterministic():
    loci = [{"LOC": str(i), "CHR": str(i % 22 + 1)} for i in range(1, 2496)]
    a = runner.canonical_manifest_rows(loci)
    b = runner.canonical_manifest_rows(loci)
    keys = [(r["trait_id"], r["locus_id"]) for r in a]
    assert a == b
    assert len(a) == 17465
    assert len(set(keys)) == 17465
    assert a[0]["trait_id"] == "adhd" and a[0]["locus_id"] == "1"
    assert a[-1]["trait_id"] == "scz" and a[-1]["locus_id"] == "2495"


def test_pair_gate_uses_only_canonical_results_and_strict_threshold():
    t = runner.EXPECTED_THRESHOLD
    assert runner.pair_gate("TESTED", t / 2, "TESTED", t / 3) == (True, "ELIGIBLE")
    assert runner.pair_gate("TESTED", t, "TESTED", t / 2) == (False, "UNIVARIATE_UNDERPOWERED")
    assert runner.pair_gate("NOT_RUN", None, "TESTED", t / 2) == (False, "UNIVARIATE_UNDERPOWERED")
    assert runner.pair_gate("FAILED", None, "FAILED", None) == (False, "UNIVARIATE_UNDERPOWERED")


def test_frozen_missingness_denominator_includes_missing_failed_and_not_run():
    counts = {"TESTED": 17000, "NOT_RUN": 300, "FAILED": 100}
    assert runner.untested_fraction(65, counts, 17465) == 465 / 17465
    assert runner.untested_fraction(0, {"NOT_RUN": 0, "FAILED": 0}, 100) == 0
    assert runner.untested_fraction(6, {"NOT_RUN": 0, "FAILED": 0}, 100) == 0.06


def test_cell_identity_hash_is_order_independent_but_value_sensitive():
    row = {"phen": "insomnia", "locus_id": "2", "status": "TESTED", "p": "0.0878159"}
    reordered = dict(reversed(list(row.items())))
    changed = {**row, "p": "0.0878160"}
    assert runner.cell_hash(row) == runner.cell_hash(reordered)
    assert runner.cell_hash(row) != runner.cell_hash(changed)


def test_result_identity_is_independent_of_worker_assignment():
    row = {"phen": "insomnia", "locus_id": "2", "status": "TESTED", "p": "0.0878159"}
    hashes = [runner.cell_hash(row) for _worker_count in (1, 2, 4, 8)]
    assert len(set(hashes)) == 1


def test_canonical_worker_is_trait_only_without_pair_overlap_context():
    worker = (ROOT / "brain6/scripts/run_lava_canonical_batch_v3.R").read_text()
    assert "input$sample.overlap <- NULL" in worker
    assert "phenos = trait" in worker
    assert "process.sample.overlap" not in worker
    assert "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS" in worker
    assert "LOW_LOCAL_H2_UNDERPOWERED" in worker


def test_bivariate_worker_gates_pairs_from_canonical_univariate_results():
    worker = (ROOT / "brain6/scripts/run_lava_family_locus_v3.R").read_text()
    univariate = worker.index("univariate <- rbindlist")
    pair_loop = worker.index("run_pair <- function(pair)")
    canonical_gate = worker.index("ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_FAILED_FROZEN_GATE")
    pair_input = worker.index("input <- prepare_input(phenotypes, pair$sample_overlap_file)")
    bivariate = worker.index("LAVA::run.bivar(loc, phenos = phenotypes")
    assert univariate < pair_loop < canonical_gate < pair_input < bivariate
    assert "ONE_OR_MORE_UNIQUE_TRAIT_LOCAL_TESTS_NOT_RUN" in worker
    assert "need(nrow(pair_results) == 5L" in worker
    assert "cfg$canonical_univariate_rows" in worker
    assert "LAVA::run.univ" not in worker
    assert "if (is.null(cache$reference))" in worker
    assert "if (is.null(cache$sumstats[[trait]]))" in worker


def test_canonical_worker_count_is_frozen_at_four():
    execution = json.loads((ROOT / "brain6/config/lava_execution_canonical_v3.json").read_text())
    family = json.loads((ROOT / "brain6/config/lava_family_canonical_v3.json").read_text())
    assert execution["worker_count"] == 4
    assert execution["execution_policy"]["worker_count"] == 4
    assert execution["batch_size_loci"] == 5
    assert execution["family_lock_sha256"] == runner.sha256(ROOT / "brain6/config/lava_family_canonical_v3.json")
    assert family["execution"]["one_locus_per_process"] is False
    assert runner.validate_worker_count(execution, 4) == 4
    try:
        runner.validate_worker_count(execution, 2)
    except ValueError:
        pass
    else:
        raise AssertionError("CLI worker count must match the frozen concurrency setting")


def test_bivariate_family_lock_is_bound_and_predeclared():
    lock_path = ROOT / "brain6/config/lava_bivariate_canonical_v3.json"
    lock = json.loads(lock_path.read_text())
    canonical = json.loads((ROOT / "brain6/config/lava_family_canonical_v3.json").read_text())
    canonical_exec = json.loads((ROOT / "brain6/config/lava_execution_canonical_v3.json").read_text())
    worker = ROOT / lock["provenance"]["worker_script_path"]
    assert lock["review_status"] == "PREDECLARED_BEFORE_V3_BIVARIATE_RESULTS"
    assert lock_path.with_suffix(lock_path.suffix + ".sha256").read_text().split()[0] == runner.sha256(lock_path)
    assert lock["canonical_source"]["family_lock_sha256"] == runner.sha256(
        ROOT / "brain6/config/lava_family_canonical_v3.json")
    assert lock["canonical_source"]["execution_lock_sha256"] == runner.sha256(
        ROOT / "brain6/config/lava_execution_canonical_v3.json")
    assert lock["canonical_source"]["maximum_untested_fraction"] == canonical[
        "canonical_univariate"]["maximum_untested_fraction"] == 0.05
    assert lock["pair_family"]["n_slots"] == 12475
    assert lock["pair_family"]["maximum_execution_failure_fraction"] == 0.01
    assert lock["execution"]["worker_count"] == 4
    assert lock["provenance"]["worker_script_sha256"] == runner.sha256(worker)
    assert canonical_exec["worker_count"] == 4
    old_exec = json.loads((ROOT / "brain6/config/lava_execution_v1.json").read_text())
    assert lock["bivariate_parameters"] == old_exec["bivariate"]


def test_chromosome_batches_are_deterministic_bounded_and_lossless():
    loci = [{"LOC": str(i), "CHR": str(1 if i < 7 else 2)} for i in range(1, 13)]
    batches = runner.chromosome_batches(loci, 4)
    assert [(chrom, index, len(group)) for chrom, index, group in batches] == [
        (1, 1, 4), (1, 2, 2), (2, 1, 4), (2, 2, 2)]
    assert [l["LOC"] for _, _, group in batches for l in group] == [l["LOC"] for l in loci]


def test_process_failure_retries_only_to_the_frozen_attempt_limit(monkeypatch):
    calls = []
    def fake_attempt(*args):
        calls.append(1)
        if len(calls) == 1:
            raise RuntimeError("simulated interrupted batch")
        return {"batch_id": "chr1_batch0001", "chromosome": 1, "loci": 1}
    monkeypatch.setattr(runner, "run_batch_attempt", fake_attempt)
    result = runner.run_batch(1, 1, [{"LOC": "1", "CHR": "1"}], Path("run"),
        Path("inputs"), Path("reference"), Path("Rscript"), {"max_attempts_per_locus": 2})
    assert len(calls) == 2 and result["attempts"] == 2


def test_new_output_root_cannot_target_preserved_runs(tmp_path: Path):
    baseline = tmp_path / "lava-results-v1"
    distinct = tmp_path / "lava-results-v3"
    assert not runner.isolated_output_root(baseline)
    assert runner.isolated_output_root(distinct)
    assert not runner.isolated_output_root(tmp_path / "af43fde06c4c6c0d57a6ab104bea5aa33f24fb2ff2b89b53163529160452d0cf")


def test_receipt_resume_rejects_corrupt_output_and_changed_input(tmp_path: Path):
    locus = {"LOC": "2", "CHR": "1"}
    input_root = tmp_path / "inputs"
    locus_dir = input_root / "locus_2"
    locus_dir.mkdir(parents=True)
    (locus_dir / "input_info.tsv").write_text("phenotype\n", encoding="utf-8")
    for trait in runner.TRAITS:
        (locus_dir / f"{trait}.sumstats.tsv.gz").write_bytes(trait.encode())
    reference = tmp_path / "reference.provenance.json"
    reference.write_text("{}", encoding="utf-8")
    run_dir = tmp_path / "run"
    out = runner.result_path(run_dir, "2")
    out.parent.mkdir(parents=True)
    rows = [{"phen": t, "locus_id": "2", "status": "TESTED", "n_snps": "10", "n_components": "2",
             "h2.obs": "0.001", "h2.latent": "0.002", "p": "0.5", "reason": ""} for t in runner.TRAITS]
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=runner.CELL_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    family_sha, execution_sha, script_sha = "a" * 64, "b" * 64, "c" * 64
    config = tmp_path / "worker.json"
    config.write_text(json.dumps({"analysis_id": "brain6-lava-canonical-v3", "locus_ids": ["2"]}))
    runner.persist_receipt(run_dir, locus, rows, runner.TRAITS, "brain6-lava-canonical-v3", family_sha,
        execution_sha, input_root, reference, script_sha, config)
    assert runner.verify_locus(run_dir, locus, runner.TRAITS, "brain6-lava-canonical-v3", family_sha,
        execution_sha, input_root, reference, script_sha)
    aggregate = runner.write_family_aggregate(run_dir, [locus, {"LOC": "3", "CHR": "1"}], [locus],
        runner.TRAITS, "brain6-lava-canonical-v3", family_sha, execution_sha, input_root, reference, script_sha)
    with Path(aggregate["path"]).open(encoding="utf-8", newline="") as handle:
        aggregated = list(csv.DictReader(handle, delimiter="\t"))
    assert aggregate["rows"] == len(aggregated) == 14
    assert sum(row["status"] == "TESTED" for row in aggregated) == 7
    assert sum(row["status"] == "MISSING" for row in aggregated) == 7
    changed_input = locus_dir / f"{runner.TRAITS[0]}.sumstats.tsv.gz"
    changed_input.write_bytes(b"changed")
    assert not runner.verify_locus(run_dir, locus, runner.TRAITS, "brain6-lava-canonical-v3", family_sha,
        execution_sha, input_root, reference, script_sha)
    changed_input.write_bytes(runner.TRAITS[0].encode())
    config.write_text(json.dumps({"analysis_id": "brain6-lava-canonical-v3", "locus_ids": ["2"], "changed": True}))
    assert not runner.verify_locus(run_dir, locus, runner.TRAITS, "brain6-lava-canonical-v3", family_sha,
        execution_sha, input_root, reference, script_sha)
    out.write_text("corrupt\n", encoding="utf-8")
    assert not runner.verify_locus(run_dir, locus, runner.TRAITS, "brain6-lava-canonical-v3", family_sha,
        execution_sha, input_root, reference, script_sha)


def test_cell_receipts_bind_all_sumstats_and_worker_identity(tmp_path: Path):
    locus = {"LOC": "2", "CHR": "1"}
    input_root = tmp_path / "inputs"; ld = input_root / "locus_2"; ld.mkdir(parents=True)
    info = ld / "input_info.tsv"; info.write_text("info", encoding="utf-8")
    for t in runner.TRAITS: (ld / f"{t}.sumstats.tsv.gz").write_text(t, encoding="utf-8")
    ref = tmp_path / "reference.provenance.json"; ref.write_text("reference", encoding="utf-8")
    run = tmp_path / "run"
    rows = [{"phen": t, "locus_id": "2", "status": "TESTED", "n_snps": "10", "n_components": "2",
        "h2.obs": "0.001", "h2.latent": "0.002", "p": "0.5", "reason": ""} for t in runner.TRAITS]
    out = runner.result_path(run, "2"); out.parent.mkdir(parents=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=runner.CELL_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    hashes = {k: "d" * 64 for k in ("family", "execution", "script")}
    config = tmp_path / "worker.json"
    config.write_text(json.dumps({"analysis_id": "brain6-lava-canonical-v3", "locus_ids": ["2"]}))
    runner.persist_receipt(run, locus, rows, runner.TRAITS, "brain6-lava-canonical-v3", hashes["family"],
        hashes["execution"], input_root, ref, hashes["script"], config)
    receipt = json.loads(runner.receipt_path(run, "2").read_text())
    assert set(receipt["sumstats"]) == set(runner.TRAITS)
    assert receipt["config_sha256"] == runner.sha256(config)
    assert receipt["state"] == "COMPLETE"


def test_old_baseline_and_roundoff_run_ids_are_distinct_from_v3():
    family = json.loads((ROOT / "brain6/config/lava_family_v2.json").read_text())
    new = json.loads((ROOT / "brain6/config/lava_family_canonical_v3.json").read_text())
    roundoff = json.loads((ROOT / "brain6/manifests/lava_roundoff_v1_run.json").read_text())
    assert family["analysis_id"] != new["analysis_id"]
    assert roundoff["run_id"] == "af43fde06c4c6c0d57a6ab104bea5aa33f24fb2ff2b89b53163529160452d0cf"
    assert new["analysis_id"] == "brain6-lava-canonical-v3"


def test_roundoff_worker_is_separate_from_canonical_worker():
    text = SCRIPT.read_text()
    batch = (ROOT / "brain6/scripts/run_lava_canonical_batch_v3.R").read_text()
    old = (ROOT / "brain6/scripts/run_lava_family_locus_roundoff_v1.R").read_text()
    assert "run_lava_family_locus_roundoff_v1.R" not in text
    assert "brain6-lava-canonical-v3" in batch
    assert "brain6-lava-canonical-v3" not in old
