#!/usr/bin/env python3
"""Dispatch one locked interpretation task from an explicit curator import record."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/interpretation_task_manifest.tsv")
    parser.add_argument("--import-dir", default="work/interpretation_imports")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest_path = root / args.manifest
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        tasks = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in tasks if row["task_id"] == args.task_id]
    if len(selected) != 1:
        fail(f"unknown or duplicate interpretation task: {args.task_id}")
    task = selected[0]
    import_path = root / args.import_dir / f"{args.task_id}.json"
    adapter_provenance = None
    if import_path.is_file():
        record = json.loads(import_path.read_text(encoding="utf-8"))
        if set(record) - {"status", "reason", "input"} or not record.get("status") or not record.get("reason"):
            fail("interpretation curator record has unexpected or missing fields")
    else:
        policy = json.loads((root / "config/interpretation_analysis_policy.json").read_text(encoding="utf-8"))
        public_pathway_sources = set(policy["public_pathway_sources"]["source_ids"])
        automatic_sources = {
            policy["promoter_mapping"]["source_id"], policy["hocomoco_v14"]["source_id"],
            policy["abc_2021"]["source_id"], policy["pchic_2016"]["source_id"],
            policy["fuma_scrna"]["source_id"], policy["catlas_adult_v4"]["source_id"],
            policy["ldsc_seg_gtex"]["source_id"],
            *policy["screen_registry_v4"]["source_ids"], *public_pathway_sources,
        }
        automatic_family = (
            task["analysis_family"] == "regulatory"
            or task["analysis_family"] == "cell_type"
            and task["source_id"] in {
                policy["fuma_scrna"]["source_id"], policy["catlas_adult_v4"]["source_id"],
                policy["ldsc_seg_gtex"]["source_id"],
            }
            or task["analysis_family"] == "pathway"
            and task["source_id"] in public_pathway_sources
        )
        if not automatic_family or task["source_id"] not in automatic_sources:
            fail(
                f"interpretation task requires an explicit curator record: {import_path.relative_to(root)}; "
                "record terminal status/reason and a real normalized input path for COMPLETED"
            )
        automatic_dir = root / "work/interpretation_automatic" / args.task_id
        automatic_result = automatic_dir / "normalized.tsv"
        adapter_provenance = automatic_dir / "adapter.provenance.json"
        adapter_script = {
            policy["hocomoco_v14"]["source_id"]: "83_run_motif_task.py",
            policy["abc_2021"]["source_id"]: "85_run_abc_task.py",
            policy["pchic_2016"]["source_id"]: "87_run_pchic_task.py",
            policy["fuma_scrna"]["source_id"]: "91_run_fuma_scrna_task.py",
            policy["catlas_adult_v4"]["source_id"]: "97_run_catlas_task.py",
            policy["ldsc_seg_gtex"]["source_id"]: "98_run_ldsc_seg_task.py",
            **{source_id: "92_run_pathway_task.py" for source_id in public_pathway_sources},
        }.get(task["source_id"], "82_run_regulatory_task.py")
        adapter_command = [
            sys.executable, str(root / "scripts" / adapter_script), args.task_id,
            "--root", str(root), "--out", str(automatic_result),
            "--provenance-out", str(adapter_provenance),
        ]
        if not args.execute:
            print("INTERPRETATION_AUTOMATIC_TASK_READY " + " ".join(adapter_command))
            return 0
        adapter_result = subprocess.run(adapter_command, cwd=root, check=False)
        if adapter_result.returncode:
            fail(f"automatic interpretation adapter failed: {args.task_id}")
        adapter = json.loads(adapter_provenance.read_text(encoding="utf-8"))
        record = {
            "status": adapter["terminal_status"], "reason": adapter["terminal_reason"],
            "input": str(automatic_result) if adapter["terminal_status"] == "COMPLETED" else None,
        }
    command = [
        sys.executable, str(root / "scripts/76_record_interpretation_task.py"), args.task_id,
        "--root", str(root), "--status", str(record["status"]), "--reason", str(record["reason"]),
    ]
    if record.get("input"):
        command.extend(["--input", str(record["input"])])
    if adapter_provenance is not None:
        command.extend(["--adapter-provenance", str(adapter_provenance)])
    if not args.execute:
        print("INTERPRETATION_TASK_READY " + " ".join(command))
        return 0
    result = subprocess.run(command, cwd=root, check=False)
    if result.returncode:
        fail(f"interpretation task recorder failed: {args.task_id}")
    print(
        f"INTERPRETATION_TASK_DISPATCHED task={args.task_id} "
        f"family={task['analysis_family']} source={task['source_id']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
