#!/usr/bin/env python3
"""Build and verify frozen MiXeR task manifests for a portable x86_64 run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file():
        fail(f"missing checksum-bound artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def build_univariate_tasks(
    panel: list[dict[str, str]], policy: dict[str, object],
) -> list[dict[str, object]]:
    tasks: list[dict[str, object]] = []
    replicates = int(policy["fit_replicates"])
    seed_offset = int(policy["replicate_seed_offset"])
    for panel_row in panel:
        trait = panel_row["trait_id"]
        prefix = f"results/mixer/univariate/{trait}"
        for replicate in range(1, replicates + 1):
            tasks.append({
                "task_id": f"MIXER_U::{trait}::rep{replicate:02d}",
                "phase": "univariate_replicate",
                "trait1": trait,
                "trait2": "",
                "replicate": replicate,
                "seed": seed_offset + replicate,
                "command": f"bash scripts/38_run_mixer_task.sh univariate {trait} {replicate}",
                "expected_outputs": ";".join((
                    f"{prefix}.fit.rep{replicate}.json",
                    f"{prefix}.test.rep{replicate}.json",
                )),
            })
        tasks.append({
            "task_id": f"MIXER_U::{trait}::combine",
            "phase": "univariate_combine",
            "trait1": trait,
            "trait2": "",
            "replicate": "",
            "seed": "",
            "command": f"bash scripts/38_run_mixer_task.sh combine-univariate {trait}",
            "expected_outputs": ";".join((
                f"{prefix}.fit.json", f"{prefix}.test.json",
                f"{prefix}.fit.summary.csv",
            )),
        })
    return tasks


def eligible_pairs(
    panel: list[dict[str, str]], univariate_rows: list[dict[str, str]],
) -> list[tuple[str, str]]:
    trait_ids = [row["trait_id"] for row in panel]
    if [row.get("trait_id") for row in univariate_rows] != trait_ids:
        fail("univariate MiXeR table differs from the ordered 45-trait panel")
    eligible = {
        row["trait_id"] for row in univariate_rows
        if row.get("bivariate_eligibility") == "ELIGIBLE"
    }
    return [
        (sleep["trait_id"], non_sleep["trait_id"])
        for sleep in panel if sleep["domain"] == "sleep" and sleep["trait_id"] in eligible
        for non_sleep in panel
        if non_sleep["domain"] != "sleep" and non_sleep["trait_id"] in eligible
    ]


def build_bivariate_tasks(
    pairs: list[tuple[str, str]], policy: dict[str, object],
) -> list[dict[str, object]]:
    tasks: list[dict[str, object]] = []
    replicates = int(policy["fit_replicates"])
    seed_offset = int(policy["replicate_seed_offset"])
    for sleep, non_sleep in pairs:
        prefix = f"results/mixer/bivariate/{sleep}_vs_{non_sleep}"
        for replicate in range(1, replicates + 1):
            tasks.append({
                "task_id": f"MIXER_B::{sleep}__{non_sleep}::rep{replicate:02d}",
                "phase": "bivariate_replicate",
                "trait1": sleep,
                "trait2": non_sleep,
                "replicate": replicate,
                "seed": seed_offset + replicate,
                "command": (
                    f"bash scripts/38_run_mixer_task.sh bivariate "
                    f"{sleep} {non_sleep} {replicate}"
                ),
                "expected_outputs": ";".join((
                    f"{prefix}.fit.rep{replicate}.json",
                    f"{prefix}.test.rep{replicate}.json",
                )),
            })
        tasks.append({
            "task_id": f"MIXER_B::{sleep}__{non_sleep}::combine",
            "phase": "bivariate_combine",
            "trait1": sleep,
            "trait2": non_sleep,
            "replicate": "",
            "seed": "",
            "command": (
                f"bash scripts/38_run_mixer_task.sh combine-bivariate {sleep} {non_sleep}"
            ),
            "expected_outputs": ";".join((
                f"{prefix}.fit.json", f"{prefix}.test.json", f"{prefix}.csv",
            )),
        })
    return tasks


def phase_paths(
    root: Path, policy: dict[str, object], phase: str,
) -> tuple[Path, Path]:
    return (
        root / str(policy[f"{phase}_task_manifest_path"]),
        root / str(policy[f"{phase}_task_lock_path"]),
    )


def expected_tasks(
    root: Path, policy: dict[str, object], panel: list[dict[str, str]], phase: str,
) -> tuple[list[dict[str, object]], str | None]:
    if phase == "univariate":
        return build_univariate_tasks(panel, policy), None
    univariate_path = root / "results/tables/mixer_univariate.tsv"
    _, univariate_rows = read_tsv(univariate_path)
    return build_bivariate_tasks(eligible_pairs(panel, univariate_rows), policy), sha256(univariate_path)


def validate_task_lock(
    root: Path, policy: dict[str, object], panel_path: Path, phase: str,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    manifest_path, lock_path = phase_paths(root, policy, phase)
    fields, rows = read_tsv(manifest_path)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"MiXeR {phase} task lock is unreadable: {exc}")
    panel_fields, panel = read_tsv(panel_path)
    if "trait_id" not in panel_fields or len(panel) != int(policy["expected_traits"]):
        fail("MiXeR task contract requires the exact locked panel")
    expected, univariate_hash = expected_tasks(root, policy, panel, phase)
    expected_text = table_text(list(policy["task_manifest_fields"]), expected)
    if (
        fields != policy["task_manifest_fields"]
        or manifest_path.read_text(encoding="utf-8") != expected_text
        or lock.get("schema_version") != "sleep-atlas-mixer-tasks.1"
        or lock.get("analysis_id") != policy["analysis_id"]
        or lock.get("phase") != phase
        or lock.get("panel_sha256") != sha256(panel_path)
        or lock.get("policy_sha256") != sha256(root / "config/mixer_analysis_policy.json")
        or lock.get("input_lock_sha256")
        != sha256(root / str(policy["input_manifest_lock_path"]))
        or lock.get("reference_provenance_sha256")
        != sha256(root / str(policy["reference_provenance_path"]))
        or lock.get("task_planner_script_sha256") != sha256(Path(__file__))
        or lock.get("task_runner_script_sha256")
        != sha256(root / "scripts/38_run_mixer_task.sh")
        or lock.get("task_manifest_sha256") != sha256(manifest_path)
        or lock.get("task_ids_in_locked_order") != [row["task_id"] for row in rows]
        or lock.get("task_count") != len(rows)
        or lock.get("univariate_table_sha256") != univariate_hash
    ):
        fail(f"MiXeR {phase} task manifest differs from its exact pre-result lock")
    return rows, lock


def execution_preflight(root: Path) -> tuple[bool, str]:
    command = [sys.executable, str(root / "scripts/35_mixer_preflight.py"), "--root", str(root)]
    result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
    return result.returncode == 0, (result.stdout + result.stderr).strip()


def write_phase(
    root: Path, policy_path: Path, policy: dict[str, object], panel_path: Path,
    phase: str, tasks: list[dict[str, object]], univariate_hash: str | None,
) -> tuple[Path, Path]:
    manifest_path, lock_path = phase_paths(root, policy, phase)
    if manifest_path.exists() or lock_path.exists():
        fail(f"immutable MiXeR {phase} task manifest or lock already exists")
    payload = table_text(list(policy["task_manifest_fields"]), tasks)
    atomic_text(manifest_path, payload)
    lock = {
        "schema_version": "sleep-atlas-mixer-tasks.1",
        "analysis_id": policy["analysis_id"],
        "phase": phase,
        "panel_sha256": sha256(panel_path),
        "policy_sha256": sha256(policy_path),
        "input_lock_sha256": sha256(root / str(policy["input_manifest_lock_path"])),
        "reference_provenance_sha256": sha256(root / str(policy["reference_provenance_path"])),
        "univariate_table_sha256": univariate_hash,
        "task_manifest": str(manifest_path.relative_to(root)),
        "task_manifest_sha256": sha256(manifest_path),
        "task_count": len(tasks),
        "task_ids_in_locked_order": [str(row["task_id"]) for row in tasks],
        "task_planner_script_sha256": sha256(Path(__file__)),
        "task_runner_script_sha256": sha256(root / "scripts/38_run_mixer_task.sh"),
    }
    atomic_json(lock_path, lock)
    return manifest_path, lock_path


def result_artifact_paths(root: Path, tasks: list[dict[str, str]]) -> list[Path]:
    relative_paths: list[str] = []
    for task in tasks:
        relative_paths.extend(path for path in task["expected_outputs"].split(";") if path)
    if len(relative_paths) != len(set(relative_paths)):
        fail("MiXeR task family has duplicate expected result artifacts")
    paths = []
    for relative_path in relative_paths:
        relative = Path(relative_path)
        if relative.is_absolute() or ".." in relative.parts:
            fail("MiXeR task family has an unsafe expected result path")
        paths.append(root / relative)
    return paths


def seal_phase_results(
    root: Path, policy: dict[str, object], panel_path: Path, phase: str,
    canonical_table_path: Path, result_rows: int,
) -> Path:
    tasks, _ = validate_task_lock(root, policy, panel_path, phase)
    provenance_path = root / str(policy[f"{phase}_result_provenance_path"])
    if provenance_path.exists():
        fail(f"immutable MiXeR {phase} result provenance already exists")
    if not canonical_table_path.is_file() or canonical_table_path.stat().st_size == 0:
        fail(f"MiXeR {phase} canonical table is absent or empty")
    artifacts = []
    for path in result_artifact_paths(root, tasks):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"MiXeR {phase} result artifact is absent or empty: {path}")
        artifacts.append({
            "path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    manifest_path, lock_path = phase_paths(root, policy, phase)
    provenance = {
        "schema_version": "sleep-atlas-mixer-results.1",
        "analysis_id": policy["analysis_id"],
        "phase": phase,
        "mixer_release": policy["mixer_release"],
        "container_digest": policy["container_amd64_manifest_digest"],
        "reference_commit": policy["mixer_reference_commit_at_lock"],
        "reference_provenance_sha256": sha256(
            root / str(policy["reference_provenance_path"])
        ),
        "input_lock_sha256": sha256(root / str(policy["input_manifest_lock_path"])),
        "task_manifest_sha256": sha256(manifest_path),
        "task_lock_sha256": sha256(lock_path),
        "task_count": len(tasks),
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
        "canonical_table": str(canonical_table_path.relative_to(root)),
        "canonical_table_bytes": canonical_table_path.stat().st_size,
        "canonical_table_sha256": sha256(canonical_table_path),
        "canonical_result_rows": result_rows,
    }
    atomic_json(provenance_path, provenance)
    return provenance_path


def validate_phase_results(
    root: Path, policy: dict[str, object], panel_path: Path, phase: str,
    canonical_table_path: Path, result_rows: int,
) -> dict[str, object]:
    tasks, _ = validate_task_lock(root, policy, panel_path, phase)
    provenance_path = root / str(policy[f"{phase}_result_provenance_path"])
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"MiXeR {phase} result provenance is unreadable: {exc}")
    manifest_path, lock_path = phase_paths(root, policy, phase)
    artifact_paths = result_artifact_paths(root, tasks)
    records = provenance.get("artifacts", [])
    if not isinstance(records, list) or len(records) != len(artifact_paths):
        fail(f"MiXeR {phase} result artifact family differs from provenance")
    for path, record in zip(artifact_paths, records):
        if (
            not isinstance(record, dict)
            or record.get("path") != str(path.relative_to(root))
            or not path.is_file() or path.stat().st_size != record.get("bytes")
            or sha256(path) != record.get("sha256")
        ):
            fail(f"MiXeR {phase} result artifact differs from provenance: {path}")
    if (
        provenance.get("schema_version") != "sleep-atlas-mixer-results.1"
        or provenance.get("analysis_id") != policy["analysis_id"]
        or provenance.get("phase") != phase
        or provenance.get("mixer_release") != policy["mixer_release"]
        or provenance.get("container_digest") != policy["container_amd64_manifest_digest"]
        or provenance.get("reference_commit") != policy["mixer_reference_commit_at_lock"]
        or provenance.get("reference_provenance_sha256")
        != sha256(root / str(policy["reference_provenance_path"]))
        or provenance.get("input_lock_sha256")
        != sha256(root / str(policy["input_manifest_lock_path"]))
        or provenance.get("task_manifest_sha256") != sha256(manifest_path)
        or provenance.get("task_lock_sha256") != sha256(lock_path)
        or provenance.get("task_count") != len(tasks)
        or provenance.get("artifact_count") != len(artifact_paths)
        or provenance.get("canonical_table") != str(canonical_table_path.relative_to(root))
        or not canonical_table_path.is_file()
        or canonical_table_path.stat().st_size != provenance.get("canonical_table_bytes")
        or sha256(canonical_table_path) != provenance.get("canonical_table_sha256")
        or provenance.get("canonical_result_rows") != result_rows
    ):
        fail(f"MiXeR {phase} result provenance differs from the locked execution contract")
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("univariate", "bivariate"))
    parser.add_argument("--root", default=".")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify-task")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy_path = root / "config/mixer_analysis_policy.json"
    panel_path = root / "config/analysis_panel.tsv"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if args.write and args.verify_task:
        fail("choose task-manifest creation or task verification, not both")
    _, panel = read_tsv(panel_path)
    if len(panel) != int(policy["expected_traits"]) or len({row["trait_id"] for row in panel}) != len(panel):
        fail("MiXeR task planning requires the exact locked 45-trait panel")
    if args.phase == "bivariate":
        validator = subprocess.run(
            [sys.executable, str(root / "scripts/40_validate_mixer.py"),
             "--root", str(root), "--univariate-only", "--quiet"],
            cwd=root, capture_output=True, text=True, check=False,
        )
        if validator.returncode:
            fail("validated complete univariate MiXeR results are required before bivariate planning")
    if args.verify_task:
        rows, _ = validate_task_lock(root, policy, panel_path, args.phase)
        selected = [row for row in rows if row["task_id"] == args.verify_task]
        if len(selected) != 1:
            fail(f"task is absent or duplicated in the locked MiXeR {args.phase} family")
        print(
            f"MIXER_TASK_VERIFIED task={args.verify_task} command={selected[0]['command']}"
        )
        return 0
    tasks, univariate_hash = expected_tasks(root, policy, panel, args.phase)
    ready, detail = execution_preflight(root)
    print(
        f"MIXER_TASK_PLAN phase={args.phase} tasks={len(tasks)} "
        f"execution_preflight={'READY' if ready else 'BLOCKED'}"
    )
    if not args.write:
        if not ready:
            print(detail)
        print("No task manifest written. Re-run with --write only on the preflight-ready x86_64 host.")
        return 0
    if not ready:
        fail("MiXeR execution preflight is blocked; no task manifest was written")
    manifest_path, lock_path = write_phase(
        root, policy_path, policy, panel_path, args.phase, tasks, univariate_hash,
    )
    validate_task_lock(root, policy, panel_path, args.phase)
    print(
        f"MIXER_TASKS_LOCKED phase={args.phase} tasks={len(tasks)} "
        f"manifest={manifest_path.relative_to(root)} lock={lock_path.relative_to(root)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
