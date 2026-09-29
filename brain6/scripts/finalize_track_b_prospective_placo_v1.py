#!/usr/bin/env python3
"""Finish the isolated Track B sensitivity when its 276 chunks are sealed."""

from __future__ import annotations

import fcntl
import json
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = Path("/Volumes/Extreme SSD/brain6-work/placo-track-b-prospective-v1")
RUNNER = ROOT / "brain6/scripts/run_track_b_prospective_placo_v1.py"
PYTHON = ROOT / "extensions/brain6/.venv/bin/python"
EXPECTED = 276
MAX_WAIT_SECONDS = 8 * 60 * 60


def completed() -> int:
    return sum((OUT / "native" / f"chunk_{i:06d}" / f"chunk_{i:06d}" / "receipt.json").is_file()
               for i in range(EXPECTED))


def main() -> None:
    contract = json.loads((OUT / "contract.json").read_text())
    if contract["analysis_id"] != "insomnia__adhd_public_input_prospective_sensitivity_v1":
        raise SystemExit("wrong analysis contract")
    lock_path = OUT / "finalization.lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.monotonic() + MAX_WAIT_SECONDS
        while completed() < EXPECTED:
            if time.monotonic() >= deadline:
                raise SystemExit(f"timed out after {MAX_WAIT_SECONDS}s with {completed()}/{EXPECTED} chunks")
            time.sleep(30)
        env = os.environ.copy()
        env["PATH"] = "/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/bin:" + env.get("PATH", "")
        with (OUT / "finalization.log").open("a") as log:
            log.write(f"all {EXPECTED} chunk receipts observed; launching receipt-verified collation\n")
            log.flush()
            result = subprocess.run(
                [str(PYTHON), str(RUNNER), "collate"],
                cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
                check=False,
            )
            if result.returncode:
                raise SystemExit(f"collation exited {result.returncode}; inspect {OUT / 'finalization.log'}")
        print(json.dumps({"status": "COLLATION_COMPLETE", "chunks": completed(),
                          "result": str(OUT / "collated/results.tsv.gz")}), flush=True)


if __name__ == "__main__":
    main()
