#!/usr/bin/env python3
"""Run one locked trait through the frozen LDSC-SEG h2-cts family."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    opener = gzip.open if path.name.endswith(".gz") else Path.open
    if opener is gzip.open:
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open(encoding="utf-8", newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def read_header(path: Path) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    if path.name.endswith(".gz"):
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open(encoding="utf-8", newline="")
    with handle:
        return handle.readline().rstrip("\n").split("\t")


def validate_result(path: Path, selected: list[dict[str, object]]) -> list[dict[str, str]]:
    fields, rows = read_tsv(path)
    expected_fields = ["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"]
    if fields != expected_fields or len(rows) != len(selected):
        fail("LDSC-SEG h2-cts output differs from its exact result family")
    expected = {str(row["source_label"]) for row in selected}
    observed: set[str] = set()
    for row in rows:
        name = row["Name"]
        if name not in expected or name in observed:
            fail(f"LDSC-SEG h2-cts output has an unexpected or duplicate tissue: {name}")
        observed.add(name)
        try:
            coefficient = float(row["Coefficient"])
            se = float(row["Coefficient_std_error"])
            p_value = float(row["Coefficient_P_value"])
        except ValueError as exc:
            fail(f"LDSC-SEG h2-cts output is nonnumeric for {name}")
            raise AssertionError from exc
        if (
            not all(math.isfinite(value) for value in (coefficient, se, p_value))
            or se <= 0 or not 0 <= p_value <= 1
        ):
            fail(f"LDSC-SEG h2-cts output is outside policy for {name}")
    if observed != expected:
        fail("LDSC-SEG h2-cts output omits a prespecified tissue")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trait_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--policy", default="config/interpretation_analysis_policy.json")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/interpretation_task_manifest.lock.json")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / args.policy
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    spec = policy["ldsc_seg_gtex"]
    selection_path = root / spec["selection_config"]
    reference_config_path = root / spec["reference_config"]
    if (
        sha256(selection_path) != spec["selection_config_sha256"]
        or sha256(reference_config_path) != spec["reference_config_sha256"]
    ):
        fail("LDSC-SEG selection or reference config differs from policy")
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_config_path.read_text(encoding="utf-8"))
    manifest_path, lock_path = root / args.manifest, root / args.manifest_lock
    fields, tasks = read_tsv(manifest_path)
    if fields != policy["task_manifest_fields"] or not lock_path.is_file():
        fail("interpretation task family is absent or has the wrong schema")
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if (
        lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("task_manifest_sha256") != sha256(manifest_path)
        or lock.get("task_ids_in_locked_order") != [row["task_id"] for row in tasks]
    ):
        fail("interpretation task family differs from its pre-result lock")
    selected_tasks = [
        row for row in tasks
        if row["analysis_family"] == "cell_type" and row["source_id"] == spec["source_id"]
        and row["trait_id"] == args.trait_id
    ]
    domains = policy["cell_types"]["required_domains"]
    if len(selected_tasks) != len(domains) or {row["domain"] for row in selected_tasks} != set(domains):
        fail("trait does not have the exact locked four-domain LDSC-SEG task family")
    reference_provenance_path = root / spec["ldscore_cache_provenance_path"]
    ldcts_path = root / spec["selected_ldcts_path"]
    if not reference_provenance_path.is_file() or not ldcts_path.is_file():
        fail("derived LDSC-SEG LD scores and selected ldcts are incomplete")
    reference_provenance = json.loads(reference_provenance_path.read_text(encoding="utf-8"))
    if (
        reference_provenance.get("reference_config_sha256") != sha256(reference_config_path)
        or reference_provenance.get("selection_config_sha256") != sha256(selection_path)
        or reference_provenance.get("selected_ldcts_sha256") != sha256(ldcts_path)
        or reference_provenance.get("derived_file_count") != 22 * 17 * 4
    ):
        fail("derived LDSC-SEG reference differs from its frozen provenance")
    for record in reference_provenance.get("derived_files", []):
        path = root / str(record.get("path", ""))
        if not path.is_file() or path.stat().st_size != record.get("bytes"):
            fail(f"derived LDSC-SEG file is absent or has the wrong size: {path}")
    sumstats = root / "data/munged" / f"{args.trait_id}.sumstats.gz"
    sumstats_fields = read_header(sumstats)
    if not {"SNP", "A1", "A2", "Z", "N"}.issubset(sumstats_fields):
        fail("LDSC-SEG input is not a canonical munged summary-statistics file")
    out_path = root / spec["trait_result_path_template"].format(trait_id=args.trait_id)
    log_path = root / spec["trait_log_path_template"].format(trait_id=args.trait_id)
    provenance_path = root / spec["trait_provenance_path_template"].format(trait_id=args.trait_id)
    if any(path.exists() for path in (out_path, log_path, provenance_path)):
        fail("LDSC-SEG trait output is immutable and already exists")
    static_root = root / reference["static_reference_dir"]
    temporary_prefix = out_path.parent / f".{args.trait_id}.ldsc-seg.tmp"
    temporary_result = Path(str(temporary_prefix) + ".cell_type_results.txt")
    temporary_log = Path(str(temporary_prefix) + ".log")
    if temporary_result.exists() or temporary_log.exists():
        fail("stale LDSC-SEG trait temporary output exists; inspect it before retrying")
    command = [
        str(root / reference["ldsc_python"]), str(root / reference["ldsc_script"]),
        "--h2-cts", str(sumstats),
        "--ref-ld-chr", str(static_root / "baselineLD_v2.2/baselineLD."),
        "--w-ld-chr", str(static_root / "1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."),
        "--ref-ld-chr-cts", str(ldcts_path), "--out", str(temporary_prefix),
    ]
    if not args.execute:
        print("LDSC_SEG_TRAIT_READY " + " ".join(command))
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
    if result.returncode or not temporary_result.is_file() or not temporary_log.is_file():
        fail(f"LDSC-SEG h2-cts failed for {args.trait_id}: {(result.stdout + result.stderr).strip()}")
    rows = validate_result(temporary_result, selection["selected_tissues"])
    temporary_result.replace(out_path)
    temporary_log.replace(log_path)
    provenance = {
        "schema_version": "sleep-atlas-ldsc-seg-trait.1",
        "analysis_id": policy["analysis_id"],
        "trait_id": args.trait_id,
        "policy_sha256": sha256(policy_path),
        "task_manifest_sha256": sha256(manifest_path),
        "task_ids": [row["task_id"] for row in selected_tasks],
        "task_input_scope_sha256": {
            row["domain"]: row["input_scope_sha256"] for row in selected_tasks
        },
        "selection_config_sha256": sha256(selection_path),
        "reference_config_sha256": sha256(reference_config_path),
        "reference_provenance_sha256": sha256(reference_provenance_path),
        "selected_ldcts_sha256": sha256(ldcts_path),
        "munged_sumstats_path": str(sumstats.relative_to(root)),
        "munged_sumstats_bytes": sumstats.stat().st_size,
        "munged_sumstats_sha256": sha256(sumstats),
        "command": command,
        "result_path": str(out_path.relative_to(root)),
        "result_rows": len(rows),
        "result_sha256": sha256(out_path),
        "log_path": str(log_path.relative_to(root)),
        "log_sha256": sha256(log_path),
        "multiple_testing_family": spec["multiple_testing_family"],
        "claim_limit": spec["claim_limit"],
        "script_sha256": sha256(Path(__file__)),
    }
    temporary_provenance = provenance_path.with_name(provenance_path.name + ".tmp")
    temporary_provenance.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    temporary_provenance.replace(provenance_path)
    print(f"LDSC_SEG_TRAIT_OK trait={args.trait_id} tissues={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
