#!/usr/bin/env python3
"""Run the locked 12 x 2,495 FI-sleep LAVA family with resumable JSON receipts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCK = ROOT / "frailty_paper/config/lava_frailty_sensitivity_v1.yaml"
PANEL = ROOT / "config/analysis_panel.tsv"
PREPARE = ROOT / "frailty_paper/scripts/43_prepare_lava_sensitivity_inputs.py"
WORKER = ROOT / "frailty_paper/scripts/44_run_lava_sensitivity_locus.R"
RUNTIME_VALIDATOR = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
R_SCRIPT = ROOT / ".r-env/bin/Rscript"
OUTPUT_ROOT = Path("/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1")
LOCUS_FILE = ROOT / "ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
TRAITS = ("accel_sleep_duration", "chronotype", "insomnia", "longsleep", "napping", "shortsleep",
          "sleep_apnea", "sleep_efficiency", "sleep_timing", "sleepdur", "sleepiness", "snoring")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_lock() -> dict:
    text = LOCK.read_text(encoding="utf-8")
    analysis = re.search(r"^analysis_id:\s*([^\s#]+)\s*$", text, re.MULTILINE)
    outputs = re.search(r"^outputs:\s*\n((?:^[ \t]+.*\n?)*)", text, re.MULTILINE)
    out = re.search(r"^[ \t]+external_run_root:\s*(.+?)\s*$", outputs.group(1), re.MULTILINE) if outputs else None
    ref_section = re.search(r"^runtime:\s*\n((?:^[ \t]+.*\n?)*)", text, re.MULTILINE)
    ref = re.search(r"^[ \t]+reference_root:\s*(.+?)\s*$", ref_section.group(1), re.MULTILINE) if ref_section else None
    if not all((analysis, out, ref)):
        raise RuntimeError("could not read required scalar values from locked YAML")
    if Path(out.group(1)) != OUTPUT_ROOT:
        raise RuntimeError("runner output root differs from analysis lock")
    return {"analysis_id": analysis.group(1), "reference_root": ref.group(1)}


def atomic_json(path: Path, payload: dict) -> None:
    tmp = path.with_name(path.name + f".{os.getpid()}.partial")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def verify_receipt(path: Path, trait: str, pair_order: int, locus_index: int) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    identity = (value.get("schema_version"), value.get("trait"), value.get("pair_order"), value.get("locus_index"))
    if identity != (1, trait, pair_order, locus_index):
        raise RuntimeError(f"receipt identity mismatch: {path}")
    if value.get("process_status") not in {"PROCESSED", "PROCESS_FAILED"}:
        raise RuntimeError(f"invalid process status in {path}")
    if value.get("univ_status") not in {"TESTED", "PARTIAL", "UNIVARIATE_FAILED", "UNIVARIATE_INVALID", "PHENOTYPE_DROPPED", "NOT_RUN"}:
        raise RuntimeError(f"invalid univariate status in {path}")
    allowed_trait_status = {"TESTED", "PROCESS_FAILED", "UNIVARIATE_FAILED", "UNIVARIATE_INVALID", "PHENOTYPE_DROPPED", "NOT_RUN"}
    if value.get("fi_status") not in allowed_trait_status or value.get("sleep_status") not in allowed_trait_status:
        raise RuntimeError(f"invalid per-trait status in {path}")
    if value.get("bivar_status") not in {"TESTED", "NOT_ELIGIBLE", "BIVARIATE_FAILED", "BIVARIATE_INVALID"}:
        raise RuntimeError(f"invalid bivariate status in {path}")
    for status_key, fields in (("fi_status", ("fi_h2_obs", "fi_p")), ("sleep_status", ("sleep_h2_obs", "sleep_p"))):
        if value[status_key] == "TESTED":
            for key in fields:
                if not isinstance(value.get(key), (float, int)):
                    raise RuntimeError(f"tested receipt lacks finite {key}: {path}")
    if value["bivar_status"] == "TESTED" and not isinstance(value.get("bivar_p"), (float, int)):
        raise RuntimeError(f"tested bivariate receipt lacks p-value: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-new", type=int, default=0,
                        help="process at most this many new receipts (0 means complete the full family)")
    args = parser.parse_args()
    if args.max_new < 0:
        parser.error("--max-new must be nonnegative")
    lock = read_lock()
    manifest_path = OUTPUT_ROOT / "input_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError("prepared-input manifest is missing; run 43_prepare_lava_sensitivity_inputs.py first")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("analysis_id") != lock["analysis_id"] or manifest.get("lock_sha256") != sha256(LOCK):
        raise RuntimeError("prepared inputs do not match the current locked analysis")
    if manifest.get("pair_count") != 12 or manifest.get("locus_count") != 2495 or manifest.get("candidate_pair_locus_slots") != 29940:
        raise RuntimeError("prepared-input family dimensions drifted")
    if set(manifest.get("pairs", {})) != set(TRAITS):
        raise RuntimeError("prepared-input pair set differs from the frozen 12-trait panel")
    subprocess.run([sys.executable, str(PREPARE), "--verify"], check=True)
    runtime = subprocess.run([str(R_SCRIPT), str(RUNTIME_VALIDATOR)], check=True, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if "status=PASS" not in runtime.stdout or "lava_commit=e729a245f7b6923967a96804fbf5246eadf2d6c6" not in runtime.stdout:
        raise RuntimeError("pinned LAVA runtime validator did not report the expected implementation")
    # The locus definition is whitespace separated and may contain repeated spaces.
    with LOCUS_FILE.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        loci = [dict(zip(header, line.split())) for line in handle if line.strip()]
    if header != ["LOC", "CHR", "START", "STOP"] or len(loci) != 2495 or len({x["LOC"] for x in loci}) != 2495:
        raise RuntimeError("frozen 2,495-locus definition failed validation")

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    lock_path = OUTPUT_ROOT / "runner.lock"
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as e:
        raise RuntimeError(f"another runner may hold {lock_path}; inspect before removing it") from e
    os.write(lock_fd, f"pid={os.getpid()} started={datetime.now(timezone.utc).isoformat()}\n".encode())
    scratch = Path(tempfile.mkdtemp(prefix="frailtylava_"))
    input_alias = scratch / "inputs"
    input_alias.symlink_to(OUTPUT_ROOT / "inputs", target_is_directory=True)
    info_root = scratch / "input_info"
    info_root.mkdir()

    processed_new = 0
    failures: list[dict] = []
    done = 0
    try:
        for pair_order, trait in enumerate(TRAITS, 1):
            pair = manifest["pairs"][trait]
            pair_root = OUTPUT_ROOT / "pairs" / trait
            receipt_dir = pair_root / "loci"
            log_dir = pair_root / "logs"
            receipt_dir.mkdir(parents=True, exist_ok=True)
            log_dir.mkdir(parents=True, exist_ok=True)
            for locus_index, locus in enumerate(loci, 1):
                receipt_path = receipt_dir / f"locus_{locus_index:04d}.json"
                if receipt_path.exists():
                    verify_receipt(receipt_path, trait, pair_order, locus_index)
                    done += 1
                    continue
                if args.max_new and processed_new >= args.max_new:
                    break
                chromosome = int(locus["CHR"])
                info_dir = info_root / trait
                info_dir.mkdir(exist_ok=True)
                info = info_dir / f"input_info_chr{chromosome:02d}.tsv"
                original_info = pair_root / f"input_info_chr{chromosome:02d}.tsv"
                with original_info.open(encoding="utf-8", newline="") as source:
                    rows = list(csv.DictReader(source, delimiter="\t"))
                if [row["phenotype"] for row in rows] != ["frailty", trait]:
                    raise RuntimeError(f"unexpected pair input.info rows: {original_info}")
                for row in rows:
                    row["filename"] = f"{row['phenotype']}/chr{chromosome:02d}.sumstats.tsv.gz"
                with info.open("w", encoding="utf-8", newline="") as output:
                    writer = csv.DictWriter(output, fieldnames=["phenotype", "cases", "controls", "prevalence", "filename"],
                                            delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
                log_path = log_dir / f"locus_{locus_index:04d}.log"
                command = [str(R_SCRIPT), str(WORKER), trait, str(pair_order), str(locus_index),
                           str(info), pair["overlap_file"], str(LOCUS_FILE), lock["reference_root"],
                           str(input_alias), str(receipt_path)]
                error = None
                for attempt in (1, 2):
                    with log_path.open("a", encoding="utf-8") as log:
                        log.write(f"\nATTEMPT {attempt} {datetime.now(timezone.utc).isoformat()}\n")
                        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
                    if receipt_path.is_file():
                        verify_receipt(receipt_path, trait, pair_order, locus_index)
                        error = None
                        break
                    error = f"worker exit {completed.returncode}; see {log_path}"
                    time.sleep(1)
                if error:
                    failure = {"trait": trait, "pair_order": pair_order, "locus_index": locus_index,
                               "LOC": locus["LOC"], "CHR": chromosome, "error": error}
                    failures.append(failure)
                    failure_path = pair_root / "launch_failures.jsonl"
                    with failure_path.open("a", encoding="utf-8") as f:
                        f.write(json.dumps(failure, sort_keys=True) + "\n")
                else:
                    done += 1
                processed_new += 1
                if processed_new % 10 == 0:
                    atomic_json(OUTPUT_ROOT / "runner_state.json", {
                        "analysis_id": lock["analysis_id"], "lock_sha256": sha256(LOCK),
                        "expected_receipts": 29940, "verified_receipts": done,
                        "launch_failures_this_invocation": len(failures), "last_trait": trait,
                        "last_locus_index": locus_index, "updated_utc": datetime.now(timezone.utc).isoformat()})
                    print(f"RUN_PROGRESS new={processed_new} verified={done} launch_failures={len(failures)}", flush=True)
            if args.max_new and processed_new >= args.max_new:
                break
        state = {"analysis_id": lock["analysis_id"], "lock_sha256": sha256(LOCK),
                 "expected_receipts": 29940, "verified_receipts": done,
                 "launch_failures_this_invocation": len(failures),
                 "updated_utc": datetime.now(timezone.utc).isoformat(),
                 "complete": done == 29940 and not failures}
        atomic_json(OUTPUT_ROOT / "runner_state.json", state)
        print(f"RUN_COMPLETE verified={done}/29940 launch_failures={len(failures)} complete={state['complete']}")
        return 0 if not failures else 2
    finally:
        os.close(lock_fd)
        lock_path.unlink(missing_ok=True)
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
