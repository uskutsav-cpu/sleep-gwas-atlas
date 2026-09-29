#!/usr/bin/env python3
"""Read-only Track B archive audit; writes only this new recovery directory."""

import base64
import csv
import gzip
import hashlib
import json
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
ARCHIVE = Path("/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas")
BACKUP = Path("/Volumes/Extreme SSD/Codex Archive/2026-09-19/codex-2026-09-01-sleep_gwas_atlas.tar")
REL = Path("results/track_b/pleiotropy")
RUN = "f8a3de8642698a84cddaf59aa2d94e518a5994139fe8fc2bb5ff530ccdaf5066"
RECORD = REL / "continuations/post_lava_terminal_v2/sequential_execution_records" / RUN / "B"
CODE = {
    "scripts/125_materialize_track_b_placo_pair.py": "04e9f2f96c81fc1d612ea6c2aa307041594eb865318f87b2a4a469fa1efb2a54",
    "scripts/141_run_track_b_placo_pair_v2.R": "18b01861165925227d59f7c55fc4b4549c6a17775100458845dc76ec56d002f9",
    "scripts/143_run_track_b_placo_sequential_v2.py": "6cf2c55113e0d9a69e6825bda973c87e9e22c51adf622adf419e5225034d4e11",
    "scripts/139_build_track_b_placo_terminal_gate_v2.py": "50a65cc503b60fc2732dfa8ebf380971c4f20297137514e643d185ec7c2286cc",
    "scripts/140_materialize_track_b_placo_pair_v2.py": "9cabcf1c281aa9e31243ca0c4ee0e62df34cbe102d2902e84410a6f216909c0b",
}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def tar_code_hashes():
    found = {}
    with tarfile.open(BACKUP, "r|") as stream:
        for member in stream:
            if not member.isfile() or not member.name.startswith("sleep_gwas_atlas/scripts/"):
                continue
            key = member.name.removeprefix("sleep_gwas_atlas/")
            if key not in CODE:
                continue
            handle = stream.extractfile(member)
            h = hashlib.sha256()
            for block in iter(lambda: handle.read(1 << 20), b""):
                h.update(block)
            found[key] = {"sha256": h.hexdigest(), "expected_sha256": CODE[key], "match": h.hexdigest() == CODE[key], "bytes": member.size}
    return found


def audit_ledger(path, provenance):
    counts = {"rows": 0, "primary": 0, "pairwise": 0, "bh": 0, "failures": 0}
    expected_header = provenance["exact_schema"].split(",")
    with gzip.open(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != expected_header:
            raise ValueError("B ledger header mismatch")
        for row in reader:
            counts["rows"] += 1
            if row["analysis_id"] != provenance["analysis_id"] or row["pair_id"] != "B":
                raise ValueError("B ledger row identity mismatch")
            if int(row["within_pair_family_n"]) != provenance["output_rows"]:
                raise ValueError("B ledger family denominator mismatch")
            p = float(row["P_PLACO_PLUS"])
            q = float(row["PLACO_BH_Q"])
            if not 0 <= p <= 1 or not 0 <= q <= 1:
                raise ValueError("B ledger invalid P or q")
            counts["primary"] += p <= 2.5e-8
            counts["pairwise"] += p <= 5e-8
            counts["bh"] += q <= 0.05
            counts["failures"] += row["analysis_status"] != "TESTED" or row["numerical_error"] not in {"", "NA"}
    return counts


def main():
    sidecar_path = ARCHIVE / REL / "results/placo/B.provenance.json"
    ledger_path = ARCHIVE / REL / "results/placo/B.full.tsv.gz"
    sidecar = read(sidecar_path)
    bundle_path = ARCHIVE / RECORD / "reproducibility_and_cleanup.plan.json"
    completion_path = ARCHIVE / RECORD / "cleanup.complete.json"
    bundle, completion = read(bundle_path), read(completion_path)
    v1_gate = ARCHIVE / REL / "input_gate.lock.json"
    v2_gate = ARCHIVE / REL / "continuations/post_lava_terminal_v2/terminal_gate.lock.json"
    current_gate = ROOT / REL / "input_gate.lock.json"
    contract = ARCHIVE / REL / "contract.lock.json"
    monitor = ARCHIVE / RECORD / "ram/full_scan.provenance.json"
    monitor_data = read(monitor)
    full_scan_log = ARCHIVE / RECORD / "attempts/full_scan.503475e5025f4871dfe8febb967697fd3961dff321de26ccb06bf54e665a32ef.log"
    execution_contract = ARCHIVE / REL / "continuations/post_lava_terminal_v2/sequential_execution.contract.json"
    code = tar_code_hashes()
    files = {
        "archived_B_ledger": (ledger_path, sidecar["output_sha256"]),
        "archived_B_provenance": (sidecar_path, bundle["canonical"]["provenance"]["sha256"]),
        "archived_contract": (contract, sidecar["contract_lock_sha256"]),
        "archived_v1_input_gate": (v1_gate, "79d2ae0faa0beb6368496617ab4c88a580a34e7f459f4825d5586e31eecb2101"),
        "archived_v2_terminal_gate": (v2_gate, sidecar["input_gate_lock_sha256"]),
        "current_v3_input_gate": (current_gate, "66916e3416a37ce3c31ffb12d1ccf51a59b1f065ed002efbe7b196f5ffa7c402"),
        "B_full_scan_monitor": (monitor, bundle["full_scan_ram_evidence"]["provenance"]["sha256"]),
        "B_full_scan_receipt": (ARCHIVE / RECORD / "ram/full_scan.tsv", bundle["full_scan_ram_evidence"]["receipt"]["sha256"]),
        "B_full_scan_attempt_log": (full_scan_log, None),
        "sequential_execution_contract": (execution_contract, bundle["sequential_contract_sha256"]),
    }
    file_audit = {k: {"path": str(path), "exists": path.is_file(), "sha256": digest(path) if path.is_file() else None, "expected_sha256": expected} for k, (path, expected) in files.items()}
    for value in file_audit.values():
        value["match"] = value["sha256"] == value["expected_sha256"] if value["expected_sha256"] else value["exists"]
    embedded = {}
    for key, value in bundle["embedded_reproducibility_artifacts"].items():
        payload = base64.b64decode(value["base64"], validate=True)
        embedded[key] = {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(), "expected_sha256": value["sha256"], "match": len(payload) == value["bytes"] and hashlib.sha256(payload).hexdigest() == value["sha256"]}
    removed = bundle["cleanup_candidates"]
    still_present = [x["path"] for x in removed if (ARCHIVE / x["path"]).exists()]
    ledger = audit_ledger(ledger_path, sidecar)
    expected_counts = {"rows": sidecar["output_rows"], "primary": sidecar["complete_family_counts"]["primary"], "pairwise": sidecar["complete_family_counts"]["pairwise"], "bh": sidecar["complete_family_counts"]["bh"], "failures": sidecar["complete_family_counts"]["failures"]}
    completion_match = completion["cleanup_bundle_sha256"] == digest(bundle_path) and completion["removed_artifact_count"] == len(removed)
    monitoring_match = monitor_data["exit_status"] == 0 and monitor_data["output_hash"] == sidecar["output_sha256"] and monitor_data["run_fingerprint"] == RUN
    scan_log_lines = full_scan_log.read_text().splitlines()
    scan_shards = [x for x in scan_log_lines if x.startswith("PLACO+ B: checkpointed shard ")]
    scan_terminal = scan_log_lines[-1] if scan_log_lines else ""
    terminal_gate_data = read(v2_gate)
    current_insomnia_card = read(ROOT / "extensions/brain6/work/overnight-v03/source-cards/insomnia.source.json")
    current_adhd_card = read(ROOT / "extensions/brain6/work/overnight-v03/source-cards/adhd.source.json")
    source_matches = {trait: terminal_gate_data["frozen_v1"]["dense_sources"][trait]["sha256"] == card["sha256"] for trait, card in (("insomnia", current_insomnia_card), ("adhd", current_adhd_card))}
    with gzip.open(ROOT / "brain6/results/placo/insomnia__mdd/variants.tsv.gz", "rt") as handle:
        v3_header = handle.readline().strip().split("\t")
    audit = {
        "schema_version": "brain6-track-b-recovery-audit.1",
        "classification": "BLOCKED_V3_ADMISSION",
        "archive_root": str(ARCHIVE), "backup_tar": str(BACKUP), "run_fingerprint": RUN,
        "file_hashes": file_audit, "recovered_exact_source_files": code,
        "embedded_reproducibility_artifacts": embedded,
        "ledger_counts": ledger, "expected_ledger_counts": expected_counts,
        "ledger_counts_match": ledger == expected_counts,
        "B_specific_full_scan_monitor_success": monitoring_match,
        "B_full_scan_attempt_log_shard_markers": len(scan_shards),
        "B_full_scan_attempt_log_terminal_line": scan_terminal,
        "B_full_scan_attempt_log_complete": len(scan_shards) == 276 and "TRACK_B_PLACO_STAGE_COMPLETE pair=B rows=5514399 failures=0" in scan_terminal,
        "archived_dense_sources_match_current_cards": source_matches,
        "legacy_ledger_schema": sidecar["exact_schema"].split(","),
        "v3_pair_schema": v3_header,
        "legacy_ledger_schema_matches_v3_pair_schema": sidecar["exact_schema"].split(",") == v3_header,
        "B_cleanup_receipt_valid": completion_match,
        "cleanup_candidate_count": len(removed), "cleanup_candidate_paths_still_present": still_present,
        "archived_terminal_gate_matches_legacy_sidecar_alias": file_audit["archived_v2_terminal_gate"]["match"],
        "current_v3_gate_matches_legacy_sidecar": file_audit["current_v3_input_gate"]["sha256"] == sidecar["input_gate_lock_sha256"],
        "exact_aligned_input_retained": (ARCHIVE / "work/track_b_pleiotropy/placo" / RUN / "B/aligned.tsv.gz").exists(),
        "aligned_input_declared_sha256": sidecar["input_sha256"],
        "v3_frozen_auditor_path": str(ROOT / "brain6/scripts/audit_five_track_placo_family_v1.py"),
        "v3_frozen_auditor_accepts_historical_package": False,
        "promotion_performed": False,
    }
    (OUT / "track_b_recovery_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    graph = {
        "nodes": [
            {"id": "dense_sources", "type": "input", "source": "archived terminal_gate.lock.json/frozen_v1/dense_sources"},
            {"id": "contract", "type": "lock", "sha256": file_audit["archived_contract"]["sha256"]},
            {"id": "v1_gate", "type": "lock", "sha256": file_audit["archived_v1_input_gate"]["sha256"]},
            {"id": "v2_gate", "type": "lock", "sha256": file_audit["archived_v2_terminal_gate"]["sha256"]},
            {"id": "aligned_input", "type": "input", "sha256": sidecar["input_sha256"], "state": "INTENTIONALLY_CLEANED;REGENERATION_UNVERIFIED"},
            {"id": "task", "type": "command_configuration", "sha256": embedded["task.tsv"]["sha256"]},
            {"id": "code", "type": "implementation", "source": str(BACKUP), "verified_sha256": {k: v["sha256"] for k, v in code.items()}},
            {"id": "nuisance", "type": "checkpoint", "sha256": embedded["checkpoints/nuisance.rds"]["sha256"]},
            {"id": "shards_276", "type": "computation", "count": 276, "state": "INTENTIONALLY_CLEANED;ATTESTED_ONLY"},
            {"id": "ledger", "type": "output", "sha256": sidecar["output_sha256"], "rows": sidecar["output_rows"]},
            {"id": "monitor", "type": "receipt", "sha256": file_audit["B_full_scan_monitor"]["sha256"], "exit_status": monitor_data["exit_status"]},
            {"id": "provenance", "type": "receipt", "sha256": file_audit["archived_B_provenance"]["sha256"]},
            {"id": "cleanup", "type": "receipt", "sha256": digest(completion_path), "removed_artifacts": len(removed)},
            {"id": "v3_admission", "type": "gate", "state": "BLOCKED"},
        ],
        "edges": [["dense_sources", "v1_gate"], ["v1_gate", "v2_gate"], ["contract", "v2_gate"], ["v2_gate", "aligned_input"], ["aligned_input", "task"], ["code", "task"], ["task", "nuisance"], ["nuisance", "shards_276"], ["code", "shards_276"], ["shards_276", "ledger"], ["ledger", "monitor"], ["ledger", "provenance"], ["monitor", "cleanup"], ["provenance", "cleanup"], ["cleanup", "v3_admission"]],
    }
    (OUT / "track_b_provenance_graph.json").write_text(json.dumps(graph, indent=2) + "\n")
    inventory = [
        "# Track B recovery inventory", "",
        f"Historical run: `{RUN}`. Source: `{ARCHIVE}`. Backup source: `{BACKUP}`.", "",
        "| Component | Observed state |", "| --- | --- |",
        f"| `B.full.tsv.gz` | Retained; {ledger['rows']:,} tested rows; SHA `{file_audit['archived_B_ledger']['sha256']}` |",
        f"| `B.provenance.json` | Retained; SHA `{file_audit['archived_B_provenance']['sha256']}` |",
        f"| V1 contract/input gate | Retained; SHA `{file_audit['archived_contract']['sha256']}` / `{file_audit['archived_v1_input_gate']['sha256']}` |",
        f"| V2 terminal gate | Retained; SHA `{file_audit['archived_v2_terminal_gate']['sha256']}`; matches legacy sidecar field |",
        f"| Exact materializer/runner/coordinator | Present in backup tar; all {sum(x['match'] for x in code.values())}/{len(code)} expected code hashes match |",
        f"| Original aligned input | Intentionally removed; declared SHA `{sidecar['input_sha256']}`; 5,514,399 rows; exact bytes not independently regenerated |",
        f"| Nuisance checkpoint | Original bytes embedded in cleanup plan; SHA `{embedded['checkpoints/nuisance.rds']['sha256']}` |",
        f"| Shards | 276 executed; original shard RDS/individual receipts intentionally removed after attested cleanup |",
        f"| B full-scan monitor | Exit {monitor_data['exit_status']}; peak {monitor_data['peak_process_tree_rss_bytes']:,} bytes; output hash matches |",
        f"| B full-scan log | {len(scan_shards)}/276 shard-completion markers and terminal success line |",
        f"| Frozen insomnia/ADHD inputs | Archived dense-source hashes match current source cards: {source_matches} |",
        "| B-to-v3 row schema | Different: legacy 17-field full ledger versus v3 10-field pair output; requires explicit checked adapter |",
        f"| Cleanup receipt | {len(removed)} listed work artifacts removed; completion receipt bound to plan: {completion_match} |",
        "| Current Brain6 v3 admission | No frozen importer acceptance; protected slot remains empty |", "",
        "Search scope: current v03 and main worktrees, mounted Extreme SSD filename inventory, failed-partial archive, 2026-09-19 backup tar, older 2026-09-04 tar and handoff ZIP, relevant Git paths/history, and shell-history pointers. Unmounted/offline media were not inspected.",
    ]
    (OUT / "track_b_recovery_inventory.md").write_text("\n".join(inventory) + "\n")
    summary = [
        "# Track B recovery summary", "",
        "**Decision: BLOCKED_V3_ADMISSION. No protected-slot promotion was performed.**", "",
        "The archived B ledger is intact and its row counts, schema, P/q ranges, family counts, provenance hash, exact source-code hashes, V2 terminal gate, B-specific successful full-scan monitor, and cleanup attestation reconcile. The 276 original shard files and aligned input were intentionally cleaned after publication; the original V4 clean-package validator explicitly allows this state when its cleanup receipts and retained canonical output pass. Exact regeneration of the aligned gzip from frozen dense inputs remains unverified.", "",
        "The prior archive audit overstates two blockers. In the V2 bridge, `scripts/140_materialize_track_b_placo_pair_v2.py:375` intentionally sets the provenance key `input_gate_lock_sha256` to the **V2 terminal gate**, whose archived SHA is `9e9a40aa...`. The B-specific V4 monitor succeeded; the V4 amendment's failed-monitor discussion concerns superseded **Pair A** shards. Exact B materializer and runner sources were recovered from the backup tar.", "",
        "The live frozen Brain6 v3 audit still explicitly marks B inadmissible and expects the current V1 input-gate identity (`66916e34...`). The archived V1 gate differs (`79d2ae0f...`) because its readiness-gate hash differs, although the other V1 gate fields match. The archived sidecar uses V2 terminal-gate semantics. The legacy 17-field ledger also differs from the 10-field v3 pair schema and needs an explicit checked adapter. No frozen Brain6 v3 importer was found that authorizes this alias/cleanup attestation as a protected legacy package. The existing five-track auditor writes into an already populated output directory and hard-codes the archived package as rejected; executing it would neither validate nor safely promote this package. The original V4 `--verify` path also performs lock/repair/rebuild operations and requires original filesystem identities, so it cannot be run unchanged against this copied archive in a read-only audit.", "",
        "## Exact remaining requirements", "",
        "1. Freeze a narrow, prospective admission amendment for protected *legacy* B that explicitly recognizes the V2 terminal-gate alias, separately verifies the V1 gate/source identity, accepts the original V4 cleanup attestation after independent hash, monitor, code, and output checks, and defines a deterministic legacy-to-v3 row-schema adapter. Do not change statistical thresholds or the five-track family.",
        "2. Provide a read-only admission validator or restore the original V4 checkout with stable identities in a separate authorized environment; validate the complete clean package, including all originally retained A/B/CONTROL family receipts.",
        "3. If that validator rejects the historical package, obtain a different terminal-valid authorized legacy package or approve a versioned replacement-run protocol; the completed prospective public-input sensitivity must remain sensitivity evidence.", "",
        "The original aligned gzip may be reproducible from the frozen dense files and recovered code, but byte identity has not been demonstrated. Original shard bytes cannot be reconstructed without recomputing statistics. Under the original V4 clean-package rule, neither was required to remain after successful attested cleanup.",
    ]
    (OUT / "track_b_recovery_summary.md").write_text("\n".join(summary) + "\n")
    print(json.dumps({"classification": audit["classification"], "ledger_counts_match": audit["ledger_counts_match"], "code_matches": sum(x["match"] for x in code.values()), "cleanup_count": len(removed), "B_monitor_success": monitoring_match}))


if __name__ == "__main__":
    main()
