#!/usr/bin/env python3
"""Fail-closed, read-only historical Track B admission preflight."""
from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/track_b_admission_v1"
CRITERIA = OUT / "criteria.lock.json"
ARCHIVE = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
BACKUP = Path("/Volumes/Extreme SSD/Codex Archive/2026-09-19/codex-2026-09-01-sleep_gwas_atlas.tar")
STAGE = Path("/Volumes/Extreme SSD/brain6-work/track-b-legacy-admission-v1")
REL = Path("results/track_b/pleiotropy")
RUN = "f8a3de8642698a84cddaf59aa2d94e518a5994139fe8fc2bb5ff530ccdaf5066"
RECORD = REL / "continuations/post_lava_terminal_v2/sequential_execution_records" / RUN / "B"

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()

def load(path: Path):
    return json.loads(path.read_text())

def table_bytes(data: bytes) -> dict[str, str]:
    rows = list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("embedded TSV must have one data row")
    return rows[0]

def readiness(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 3:
        raise ValueError("readiness gate row count changed")
    return {r["component_id"]: r for r in rows}

def main() -> None:
    criteria = load(CRITERIA)
    expected = criteria["expected"]
    results: dict[str, dict] = {}
    def record(name: str, predicate: bool, evidence=None):
        results[name] = {"pass": bool(predicate), "evidence": evidence}
    def same(name: str, path: Path, expected_hash: str):
        if not path.is_file():
            record(name, False, {"path": str(path), "missing": True})
            return False
        actual = sha(path)
        record(name, actual == expected_hash,
               {"path": str(path), "sha256": actual, "expected_sha256": expected_hash})
        return actual == expected_hash

    auth = load(OUT / "authorization.json")
    proposal = Path(auth["proposal_path"])
    record("authorization_and_proposal_identity",
           auth["status"] == "USER_AUTHORIZED_CONDITIONAL_PROMOTION"
           and auth["proposal_sha256"] == criteria["proposal_sha256"] == sha(proposal)
           and criteria["authorization_sha256"] == sha(OUT / "authorization.json")
           and criteria["historical_run_fingerprint"] == RUN
           and criteria["promotion_rule"].startswith("ALL_REQUIRED_CHECKS_TRUE"),
           {"proposal_sha256": sha(proposal), "criteria_sha256": sha(CRITERIA)})

    files = {
        "legacy_ledger_hash": (ARCHIVE / REL / "results/placo/B.full.tsv.gz", expected["ledger_sha256"]),
        "legacy_provenance_hash": (ARCHIVE / REL / "results/placo/B.provenance.json", expected["provenance_sha256"]),
        "historical_contract_hash": (ARCHIVE / REL / "contract.lock.json", expected["contract_sha256"]),
        "historical_v1_gate_hash": (ARCHIVE / REL / "input_gate.lock.json", expected["historical_v1_gate_sha256"]),
        "historical_v2_terminal_gate_hash": (ARCHIVE / REL / "continuations/post_lava_terminal_v2/terminal_gate.lock.json", expected["historical_v2_terminal_gate_sha256"]),
        "current_v1_gate_hash": (ROOT / REL / "input_gate.lock.json", expected["current_v1_gate_sha256"]),
        "cleanup_bundle_hash": (ARCHIVE / RECORD / "reproducibility_and_cleanup.plan.json", expected["cleanup_bundle_sha256"]),
        "cleanup_completion_hash": (ARCHIVE / RECORD / "cleanup.complete.json", expected["cleanup_completion_sha256"]),
        "B_full_scan_monitor_hash": (ARCHIVE / RECORD / "ram/full_scan.provenance.json", expected["full_scan_monitor_sha256"]),
        "B_full_scan_receipt_hash": (ARCHIVE / RECORD / "ram/full_scan.tsv", expected["full_scan_receipt_sha256"]),
        "sequential_contract_hash": (ARCHIVE / REL / "continuations/post_lava_terminal_v2/sequential_execution.contract.json", expected["sequential_contract_sha256"]),
        "v3_family_lock_hash": (ROOT / "brain6/config/placo_family_v3/family_lock.json", expected["family_lock_sha256"]),
    }
    file_results = [same(name, path, value) for name, (path, value) in files.items()]
    all_files_match = all(file_results)
    record("all_exact_file_hashes", all_files_match, list(files))
    if not all_files_match:
        finish(criteria, results)
        return

    sidecar = load(files["legacy_provenance_hash"][0])
    family = load(files["v3_family_lock_hash"][0])
    old_gate = load(files["historical_v1_gate_hash"][0])
    new_gate = load(files["current_v1_gate_hash"][0])
    terminal = load(files["historical_v2_terminal_gate_hash"][0])
    bundle = load(files["cleanup_bundle_hash"][0])
    completion = load(files["cleanup_completion_hash"][0])
    monitor = load(files["B_full_scan_monitor_hash"][0])
    record("historical_v1_v2_gate_lineage",
           sidecar.get("contract_lock_sha256") == expected["contract_sha256"]
           and sidecar.get("input_gate_lock_sha256") == expected["historical_v2_terminal_gate_sha256"]
           and bundle.get("terminal_gate_lock_sha256") == expected["historical_v2_terminal_gate_sha256"]
           and old_gate.get("contract_lock_sha256") == expected["contract_sha256"]
           and new_gate.get("contract_lock_sha256") == expected["contract_sha256"]
           and terminal.get("frozen_v1", {}).get("input_gate_lock", {}).get("sha256") == expected["historical_v1_gate_sha256"],
           {"legacy_sidecar_input_gate_key": sidecar.get("input_gate_lock_sha256"),
            "historical_v1_gate": expected["historical_v1_gate_sha256"]})

    old_ready = ARCHIVE / REL / "readiness_gate.tsv"
    new_ready = ROOT / REL / "readiness_gate.tsv"
    a, b = readiness(old_ready), readiness(new_ready)
    changed = []
    for component in sorted(a):
        for field in a[component]:
            if a[component][field] != b[component][field]:
                changed.append((component, field, a[component][field], b[component][field]))
    gate_copy = {k:v for k,v in new_gate.items() if k != "readiness_gate_sha256"}
    old_copy = {k:v for k,v in old_gate.items() if k != "readiness_gate_sha256"}
    record("exact_readiness_gate_delta",
           old_copy == gate_copy and old_gate["readiness_gate_sha256"] == sha(old_ready)
           and new_gate["readiness_gate_sha256"] == sha(new_ready)
           and len(changed) == 2
           and {x[0] for x in changed} == {"PLACO_PLUS", "METHOD_UNION_08_09"}
           and all(x[1] == "blockers" and x[3] == x[2] + ";ACTIVE_LAVA_DOWNLOAD_BLOCKS_LD_MATERIALIZATION" for x in changed),
           {"changed_cells": [{"component":x[0], "field":x[1]} for x in changed]})

    source_cards = {trait: load(ROOT / f"extensions/brain6/work/overnight-v03/source-cards/{trait}.source.json") for trait in ("insomnia", "adhd")}
    frozen_sources = terminal["frozen_v1"]["dense_sources"]
    dense_lock_match = all(old_gate["live_dense_inputs"][trait]["sha256"] == new_gate["live_dense_inputs"][trait]["sha256"]
                           == frozen_sources[trait]["sha256"] for trait in old_gate["live_dense_inputs"])
    record("five_dense_source_lock_identities", dense_lock_match,
           {"traits": sorted(old_gate["live_dense_inputs"])})
    dense_bytes = {trait: same(f"{trait}_dense_bytes", Path(card["path"]), expected[f"{trait}_dense_sha256"])
                   for trait, card in source_cards.items()}
    record("insomnia_adhd_dense_bytes",
           all(dense_bytes.values())
           and all(source_cards[t]["sha256"] == expected[f"{t}_dense_sha256"]
                   == frozen_sources[t]["sha256"] for t in source_cards),
           dense_bytes)

    stage_code = {}
    for name, expected_hash in expected["source_code_sha256"].items():
        stage_code[name] = same("staged_code_" + Path(name).name, STAGE / name, expected_hash)
    found = {}
    with tarfile.open(BACKUP, "r|") as stream:
        for member in stream:
            if not member.isfile():
                continue
            name = member.name.removeprefix("sleep_gwas_atlas/")
            if name not in expected["source_code_sha256"]:
                continue
            f = stream.extractfile(member)
            h = hashlib.sha256()
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
            found[name] = h.hexdigest()
    record("recovered_code_hashes_and_tar_members",
           all(stage_code.values()) and found == expected["source_code_sha256"],
           {"tar_member_sha256": found})

    placo_path = Path(family["placo_source_path"])
    placo_ok = same("official_placo_source_bytes", placo_path, expected["placo_source_sha256"])
    historical_runner = (STAGE / "scripts/141_run_track_b_placo_pair_v2.R").read_text()
    modern_runner = (ROOT / "extensions/brain6/brain6/resources/placo.R").read_text()
    policy_path = ROOT / "config/track_b_pleiotropy_policy.json"
    policy = load(policy_path)
    record("official_placo_source_hash_and_method",
           placo_ok
           and family["placo_source_sha256"] == policy["placo_plus"]["source_sha256"] == expected["placo_source_sha256"]
           and sidecar["software_sha256_or_commit"] == expected["placo_source_sha256"]
           and all(x in historical_runner and x in modern_runner for x in ("var.placo", "cor.pearson", "placo.plus"))
           and sidecar["nuisance_estimation"]["method"] == "PINNED_PLACO_PLUS_var.placo_AND_cor.pearson"
           and sidecar["nuisance_estimation"]["marginal_p_threshold"] == criteria["frozen_statistics"]["marginal_p_threshold"]
           and policy["placo_plus"]["z_squared_maximum"] == criteria["frozen_statistics"]["z_squared_maximum"]
           and policy["placo_plus"]["absolute_tolerance"] == criteria["frozen_statistics"]["historical_absolute_tolerance"]
           and family["absolute_tolerance"] == criteria["frozen_statistics"]["v3_absolute_tolerance"],
           {"historical_abs_tol":policy["placo_plus"]["absolute_tolerance"],
            "v3_abs_tol":family["absolute_tolerance"],
            "same_official_source_sha256":expected["placo_source_sha256"]})

    embedded = bundle.get("embedded_reproducibility_artifacts", {})
    decoded = {}
    embedded_ok = True
    for name, info in embedded.items():
        try:
            data = base64.b64decode(info["base64"], validate=True)
            good = len(data) == info["bytes"] and hashlib.sha256(data).hexdigest() == info["sha256"]
            decoded[name] = data
            embedded_ok &= good
        except (ValueError, KeyError):
            embedded_ok = False
    required_embedded = {"task.tsv", "run.summary.tsv", "checkpoints/nuisance.rds",
                         "checkpoints/nuisance.sha256.tsv", "materialization.provenance.json",
                         "benchmark.raw.tsv", "benchmark.tsv"}
    embedded_ok &= set(decoded) == required_embedded
    task = table_bytes(decoded["task.tsv"]) if embedded_ok else {}
    summary = table_bytes(decoded["run.summary.tsv"]) if embedded_ok else {}
    nuisance_lock = table_bytes(decoded["checkpoints/nuisance.sha256.tsv"]) if embedded_ok else {}
    nuisance = sidecar["nuisance_estimation"]
    task_ok = embedded_ok and (
        task["run_fingerprint"] == summary["run_fingerprint"] == nuisance_lock["run_fingerprint"] == RUN
        and task["pair_id"] == summary["pair_id"] == "B"
        and task["trait1"] == sidecar["trait1"] == "insomnia"
        and task["trait2"] == sidecar["trait2"] == "adhd"
        and task["aligned_input_sha256"] == summary["aligned_sha256"] == nuisance_lock["input_sha256"] == sidecar["input_sha256"]
        and int(task["aligned_input_rows"]) == int(summary["aligned_rows"]) == expected["rows"]
        and task["input_gate_lock_sha256"] == expected["historical_v2_terminal_gate_sha256"]
        and task["contract_lock_sha256"] == expected["contract_sha256"]
        and task["placo_source_sha256"] == expected["placo_source_sha256"]
        and task["materializer_sha256"] == summary["materializer_sha256"] == expected["source_code_sha256"]["scripts/125_materialize_track_b_placo_pair.py"]
        and task["runner_sha256"] == summary["runner_sha256"] == expected["source_code_sha256"]["scripts/141_run_track_b_placo_pair_v2.R"]
        and summary["task_sha256"] == nuisance_lock["task_sha256"] == sidecar["task_sha256"] == embedded["task.tsv"]["sha256"]
        and summary["nuisance_checkpoint_sha256"] == nuisance_lock["nuisance_rds_sha256"] == nuisance["checkpoint_sha256"] == embedded["checkpoints/nuisance.rds"]["sha256"]
        and float(summary["VarZ1"]) == nuisance["VarZ1"]
        and float(summary["VarZ2"]) == nuisance["VarZ2"]
        and float(summary["CorZ"]) == nuisance["CorZ"]
        and int(summary["nuisance_null_rows_variance"]) == nuisance["variance_null_rows"]
        and int(summary["nuisance_null_rows_correlation"]) == nuisance["correlation_null_rows"]
        and float(task["marginal_p_threshold"]) == criteria["frozen_statistics"]["marginal_p_threshold"]
        and float(task["z_squared_maximum"]) == criteria["frozen_statistics"]["z_squared_maximum"]
        and float(task["absolute_tolerance"]) == criteria["frozen_statistics"]["historical_absolute_tolerance"]
        and int(summary["shard_count"]) == expected["shards"]
        and summary["terminal_status"] == sidecar["terminal_result_state"] == "COMPLETE_WITH_HITS"
    )
    record("run_task_nuisance_checkpoint_identity", task_ok,
           {"embedded_artifact_count":len(decoded), "task_sha256":embedded.get("task.tsv",{}).get("sha256"),
            "nuisance_sha256":nuisance.get("checkpoint_sha256")})

    scan_log_matches = []
    attempt_dir = ARCHIVE / RECORD / "attempts"
    for path in attempt_dir.glob("full_scan.*.log"):
        lines = path.read_text().splitlines()
        markers = sum(x.startswith("PLACO+ B: checkpointed shard ") for x in lines)
        terminal_line = lines[-1] if lines else ""
        if markers == expected["shards"] and f"TRACK_B_PLACO_STAGE_COMPLETE pair=B rows={expected['rows']} failures=0" in terminal_line:
            scan_log_matches.append(str(path))
    receipt_path = files["B_full_scan_receipt_hash"][0]
    with receipt_path.open(newline="") as f:
        ram_rows = list(csv.DictReader(f, delimiter="\t"))
    record("B_full_scan_monitor_exit_output_and_276_shards",
           monitor.get("exit_status") == 0
           and monitor.get("output_hash") == expected["ledger_sha256"]
           and monitor.get("run_fingerprint") == RUN
           and monitor.get("phase") == "full_scan"
           and len(scan_log_matches) == 1 and len(ram_rows) == 1
           and bundle["full_scan_ram_evidence"]["output_hash"] == expected["ledger_sha256"]
           and bundle["full_scan_ram_evidence"]["provenance"]["sha256"] == expected["full_scan_monitor_sha256"]
           and bundle["full_scan_ram_evidence"]["receipt"]["sha256"] == expected["full_scan_receipt_sha256"],
           {"successful_full_scan_log":scan_log_matches, "exit_status":monitor.get("exit_status"),
            "peak_process_tree_rss_bytes":monitor.get("peak_process_tree_rss_bytes")})

    retained = [bundle["canonical"]["ledger"], bundle["canonical"]["provenance"]]
    retained += bundle.get("ram_measurement_provenance", []) + bundle.get("ram_receipts", [])
    retained_status = {}
    for item in retained:
        path = ARCHIVE / item["path"]
        retained_status[item["path"]] = path.is_file() and path.stat().st_size == item["bytes"] and sha(path) == item["sha256"]
    cleanup_candidates = bundle.get("cleanup_candidates", [])
    removed_present = [x["path"] for x in cleanup_candidates if (ARCHIVE / x["path"]).exists()]
    deep = bundle["deep_validation"]
    deep_required = ("all_numerical_failures_retained_with_p_one", "bh_validated_once_over_complete_pair_family",
                     "complete_aligned_family", "fresh_global_nuisance_in_matched_benchmark",
                     "global_nuisance_all_valid_variants", "successful_full_scan_ram_evidence")
    record("V4_cleanup_bundle_and_all_retained_receipts",
           all(retained_status.values())
           and completion.get("cleanup_bundle_sha256") == expected["cleanup_bundle_sha256"]
           and completion.get("removed_artifact_count") == expected["cleanup_candidates"]
           and bundle.get("run_fingerprint") == RUN and bundle.get("pair_id") == "B"
           and bundle.get("coordinator_sha256") == expected["source_code_sha256"]["scripts/143_run_track_b_placo_sequential_v2.py"]
           and bundle.get("source_or_canonical_deletion") is False
           and all(deep.get(x) is True for x in deep_required)
           and deep["counts"] == sidecar["complete_family_counts"]
           and embedded_ok,
           {"retained_file_hashes_pass":retained_status, "embedded_count":len(decoded),
            "deep_validation":{x:deep.get(x) for x in deep_required}})
    record("all_deleted_work_artifacts_absent",
           len(cleanup_candidates) == expected["cleanup_candidates"] and not removed_present,
           {"candidate_count":len(cleanup_candidates), "unexpected_present_paths":removed_present[:10]})

    ledger_path = files["legacy_ledger_hash"][0]
    record("legacy_sidecar_result_identity",
           sidecar.get("analysis_id") == "track-b-v1.0-pleiotropy"
           and sidecar.get("pair_id") == "B"
           and sidecar.get("run_fingerprint") == RUN
           and sidecar.get("output_sha256") == expected["ledger_sha256"]
           and sidecar.get("output_rows") == expected["rows"]
           and sidecar.get("output_bytes") == ledger_path.stat().st_size
           and sidecar.get("qc_status") == "PASS"
           and sidecar.get("terminal_result_state") == "COMPLETE_WITH_HITS"
           and sidecar.get("within_pair_bh_family_n") == expected["rows"]
           and sidecar.get("policy_sha256") == sha(policy_path)
           and sidecar.get("materialization_provenance_sha256") == embedded["materialization.provenance.json"]["sha256"]
           and sidecar.get("runner_sha256") == task.get("runner_sha256")
           and sidecar.get("materializer_sha256") == task.get("materializer_sha256"),
           {"rows":sidecar.get("output_rows"), "qc_status":sidecar.get("qc_status")})

    v3_qc = load(ROOT / "brain6/results/placo/placo_v3_pair_qc_validation.json")
    four = {}
    for pair, digest in v3_qc["published_output_sha256"].items():
        path = ROOT / "brain6/results/placo" / pair / "variants.tsv.gz"
        four[pair] = path.is_file() and sha(path) == digest
    record("four_v3_pair_receipts_unchanged",
           len(four) == 4 and all(four.values())
           and v3_qc["family_lock_sha256"] == expected["family_lock_sha256"]
           and v3_qc["full_family_complete"] is False,
           four)
    finish(criteria, results)

def finish(criteria: dict, results: dict) -> None:
    required_preflight = [x for x in criteria["required_checks"]
                          if x not in {"ledger_schema_row_numeric_allele_coordinate_uniqueness",
                                       "full_5514399_BH_q_independent_check",
                                       "deterministic_17_to_10_field_adapter", "adapted_result_validation"}]
    mapping = {
        "authorization_and_proposal_identity":"authorization_and_proposal_identity",
        "all_exact_file_hashes":"all_exact_file_hashes",
        "historical_v1_v2_gate_lineage":"historical_v1_v2_gate_lineage",
        "exact_readiness_gate_delta":"exact_readiness_gate_delta",
        "five_dense_source_lock_identities":"five_dense_source_lock_identities",
        "insomnia_adhd_dense_bytes":"insomnia_adhd_dense_bytes",
        "recovered_code_hashes_and_tar_members":"recovered_code_hashes_and_tar_members",
        "official_placo_source_hash_and_method":"official_placo_source_hash_and_method",
        "run_task_nuisance_checkpoint_identity":"run_task_nuisance_checkpoint_identity",
        "B_full_scan_monitor_exit_output_and_276_shards":"B_full_scan_monitor_exit_output_and_276_shards",
        "V4_cleanup_bundle_and_all_retained_receipts":"V4_cleanup_bundle_and_all_retained_receipts",
        "all_deleted_work_artifacts_absent":"all_deleted_work_artifacts_absent",
        "four_v3_pair_receipts_unchanged":"four_v3_pair_receipts_unchanged",
    }
    passed = (all(results.get(mapping[x], {}).get("pass") is True for x in required_preflight)
              and all(item.get("pass") is True for item in results.values()))
    payload = {"analysis_id":criteria["analysis_id"], "status":"PASS" if passed else "FAIL",
               "criteria_sha256":sha(CRITERIA),"validator_sha256":sha(Path(__file__)),
               "required_preflight_checks":required_preflight,"checks":results,
               "failed_checks":[x for x in required_preflight if not results.get(mapping[x],{}).get("pass")]
                               + [x for x, item in results.items()
                                  if not item.get("pass") and x not in required_preflight]}
    path = OUT / "preflight.json"
    if path.exists() and path.read_text() != json.dumps(payload,indent=2,sort_keys=True)+"\n":
        raise RuntimeError("refusing to overwrite a different preflight result")
    if not path.exists():
        path.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":payload["status"],"failed_checks":payload["failed_checks"]}))
    if not passed:
        raise SystemExit(2)

if __name__ == "__main__":
    main()
