#!/usr/bin/env python3
"""Resume a source-specific Brain6 rescue through frozen, receipt-gated stages.

No stage command is supplied by this repository before the external source has a
reviewed statistical mapping. Stage adapters must be versioned, hashed Python
scripts in brain6/scripts; commands are argv arrays, never shell strings.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from brain6_confirmatory_restart_v1 import ROOT, preflight, sha


STAGES = (
    "harmonization",
    "frozen_88_locus_pilot",
    "full_family_projection",
    "high_impact_trait_repairs",
    "full_lava_family",
    "frozen_qc_gate",
    "candidate_promotion",
    "confirmatory_downstream",
    "regenerated_final_tables",
)
SELECTION_SHA = "9182ebca331a7362d96359cc41e130571ce274f7869477864ce7ca0a944a01e7"


def atomic_new_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    if path.exists():
        tmp.unlink()
        raise FileExistsError(path)
    os.replace(tmp, path)


def require_stage_receipt(stage: str, receipt: Path, analysis_id: str, admission_sha: str) -> dict:
    if not receipt.is_file():
        raise ValueError(f"{stage}: stage receipt missing")
    data = json.loads(receipt.read_text())
    if data.get("status") != "PASS" or data.get("stage") != stage or data.get("analysis_id") != analysis_id:
        raise ValueError(f"{stage}: receipt status/stage/analysis identity invalid")
    if data.get("source_admission_sha256") != admission_sha:
        raise ValueError(f"{stage}: receipt is not bound to frozen source admission")
    if stage == "frozen_88_locus_pilot":
        if data.get("locus_selection_sha256") != SELECTION_SHA or data.get("selected_loci") != 88 or data.get(
            "worker_count") != 4:
            raise ValueError("88-locus pilot did not use frozen selection")
        if data.get("candidate_tested", -1) - data.get("canonical_tested", -1) < 9:
            raise ValueError("88-locus pilot did not meet frozen nine-locus net-gain rule")
    if stage in ("full_family_projection", "high_impact_trait_repairs"):
        if data.get("projected_not_run_upper_bound", 17466) > 873 or data.get("source_valid") is not True:
            raise ValueError(f"{stage}: source-valid full-family route to ≤873 NOT_RUN not shown")
    if stage == "full_lava_family":
        if data.get("planned_cells") != 17465 or data.get("observed_cells") != 17465 or data.get("worker_count") != 4:
            raise ValueError("full LAVA family did not account for all 17,465 intended cells")
    if stage == "frozen_qc_gate":
        if data.get("planned_cells") != 17465 or data.get("not_run_cells", 17466) > 873 or data.get(
            "frozen_gate_pass") is not True:
            raise ValueError("frozen full-family QC gate failed")
    if stage in ("candidate_promotion", "regenerated_final_tables"):
        if data.get("pair_specific_candidate_rows") != 25 or data.get("geographic_regions") != 20:
            raise ValueError(f"{stage}: protected candidate/region counts changed")
    return data


def stage_binding(manifest: dict, stage: str, output: Path) -> tuple[Path, Path, list[str]]:
    spec = manifest.get("workflow_stages", {}).get(stage)
    if not isinstance(spec, dict):
        raise ValueError(f"{stage}: versioned command/receipt binding missing")
    script = Path(str(spec.get("script", "")))
    if script.is_absolute():
        raise ValueError(f"{stage}: script must be repository-relative")
    script = (ROOT / script).resolve()
    if not script.is_relative_to(ROOT / "brain6/scripts") or not script.is_file() or script.suffix != ".py":
        raise ValueError(f"{stage}: script must be a versioned brain6/scripts Python file")
    if sha(script) != spec.get("script_sha256"):
        raise ValueError(f"{stage}: script SHA-256 mismatch")
    relative_receipt = Path(str(spec.get("receipt", "")))
    receipt = (output / relative_receipt).resolve()
    if relative_receipt.is_absolute() or not receipt.is_relative_to(output) or receipt == output:
        raise ValueError(f"{stage}: receipt must be inside the new isolated output root")
    args = spec.get("args")
    if not isinstance(args, list) or any(not isinstance(x, str) for x in args):
        raise ValueError(f"{stage}: args must be a string argv array")
    replacements = {"$OUTPUT_ROOT": str(output), "$ANALYSIS_ID": str(manifest["analysis_id"])}
    argv = [sys.executable, str(script), *(replacements.get(x, x) for x in args)]
    return script, receipt, argv


def run(manifest_path: Path, execute: bool) -> int:
    manifest = json.loads(manifest_path.read_text())
    output = Path(str(manifest.get("new_output_root", ""))).expanduser().resolve()
    errors, freeze = preflight(manifest, manifest_path, allow_existing_output=True)
    if errors:
        print(json.dumps({"status": "BLOCKED_SOURCE_ADMISSION", "errors": errors}, indent=2))
        return 2
    admission = output / "source_admission.freeze.json"
    if admission.exists():
        if json.loads(admission.read_text()) != freeze:
            raise ValueError("Existing admission receipt differs from current source/manifest/anchors")
    elif execute:
        if output.exists():
            raise ValueError("Output root exists without a matching admission receipt")
        output.mkdir(parents=True, exist_ok=False)
        atomic_new_json(admission, freeze)
    admission_sha = sha(admission) if admission.exists() else "PENDING_FREEZE"
    progress = []
    for stage in STAGES:
        try:
            _script, receipt, argv = stage_binding(manifest, stage, output)
        except ValueError as exc:
            progress.append({"stage": stage, "status": "BLOCKED", "reason": str(exc)})
            break
        checkpoint = output / "workflow_checkpoints" / f"{stage}.json"
        if checkpoint.exists():
            saved = json.loads(checkpoint.read_text())
            if not receipt.is_file() or sha(receipt) != saved.get("receipt_sha256"):
                raise ValueError(f"{stage}: checkpoint/receipt mismatch")
            require_stage_receipt(stage, receipt, manifest["analysis_id"], admission_sha)
            progress.append({"stage": stage, "status": "VERIFIED_COMPLETE"})
            continue
        if receipt.exists():
            if not admission.exists():
                raise ValueError(f"{stage}: cannot adopt receipt before source admission")
            require_stage_receipt(stage, receipt, manifest["analysis_id"], admission_sha)
            if execute:
                atomic_new_json(checkpoint, {"stage": stage, "receipt_sha256": sha(receipt),
                                             "script_sha256": sha(Path(argv[1])), "argv": argv})
            progress.append({"stage": stage, "status": "VERIFIED_RECEIPT"})
            continue
        if not execute:
            progress.append({"stage": stage, "status": "READY_AFTER_PRIOR_GATES", "argv": argv})
            continue
        if not admission.exists():
            raise ValueError("Source admission not frozen")
        completed = subprocess.run(argv, cwd=ROOT, check=False)
        if completed.returncode != 0:
            progress.append({"stage": stage, "status": "FAILED_COMMAND", "exit_code": completed.returncode})
            break
        require_stage_receipt(stage, receipt, manifest["analysis_id"], admission_sha)
        atomic_new_json(checkpoint, {"stage": stage, "receipt_sha256": sha(receipt),
                                     "script_sha256": sha(Path(argv[1])), "argv": argv})
        progress.append({"stage": stage, "status": "VERIFIED_COMPLETE"})
    done = len(progress) == len(STAGES) and all(p["status"] in ("VERIFIED_COMPLETE", "VERIFIED_RECEIPT") for p in progress)
    print(json.dumps({"status": "COMPLETE" if done else "BLOCKED_OR_PENDING", "analysis_id": manifest["analysis_id"],
                      "source_admission": "FROZEN" if admission.exists() else "READY_TO_FREEZE", "stages": progress}, indent=2))
    return 0 if done or (not execute and all(p["status"] != "BLOCKED" for p in progress)) else 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--execute", action="store_true", help="run only source-admitted, hashed, receipt-gated stages")
    args = parser.parse_args()
    return run(args.manifest, args.execute)


if __name__ == "__main__":
    raise SystemExit(main())
