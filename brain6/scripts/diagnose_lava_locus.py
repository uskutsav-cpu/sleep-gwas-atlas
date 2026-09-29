#!/usr/bin/env python3
"""Replay a completed LAVA locus using its receipt-bound worker configuration."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions" / "brain6"))
sys.path.insert(0, str(ROOT / "brain6" / "scripts"))
from run_lava_family import EXPECTED_PAIRS, verify_unit  # noqa: E402
from brain6.io import read_json  # noqa: E402


def diagnose(locus_dir: Path, pair_id: str | None = None, rscript: str = "Rscript") -> int:
    locus_dir = locus_dir.resolve()
    receipt_path = locus_dir / "receipt.json"
    receipt = read_json(receipt_path)
    locus_id = str(receipt.get("locus_id", ""))
    run_id = str(receipt.get("run_id", ""))
    if not locus_id or not run_id:
        raise ValueError("receipt must contain locus_id and run_id")
    if verify_unit(locus_dir, locus_id, set(EXPECTED_PAIRS), run_id) is None:
        raise ValueError(f"no complete, verified receipt in {locus_dir}")
    config_path = locus_dir / "worker_config.json"
    cfg = read_json(config_path)
    configured = {p["pair_id"] for p in cfg.get("pairs", [])}
    if configured != set(EXPECTED_PAIRS):
        raise ValueError("worker configuration pair inventory does not match the frozen family")
    if pair_id is not None and pair_id not in configured:
        raise ValueError(f"pair is not in the frozen family: {pair_id}")
    # Pass a transient, read-only selection file outside the production attempt.
    selected = dict(cfg)
    selected["pairs"] = [p for p in cfg["pairs"] if pair_id is None or p["pair_id"] == pair_id]
    import tempfile
    with tempfile.TemporaryDirectory(prefix="brain6-lava-diagnostic-") as tmp:
        diagnostic_config = Path(tmp) / "config.json"
        diagnostic_config.write_text(json.dumps(selected), encoding="utf-8")
        script = ROOT / "brain6/scripts/diagnose_lava_locus.R"
        result = subprocess.run([rscript, str(script), str(diagnostic_config)],
                                cwd=ROOT, check=False)
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("locus_dir", type=Path, help="completed locus attempt directory")
    parser.add_argument("--pair", dest="pair_id", help="diagnose one frozen pair (default: all five)")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    args = parser.parse_args()
    try:
        return diagnose(args.locus_dir, args.pair_id, args.rscript)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
