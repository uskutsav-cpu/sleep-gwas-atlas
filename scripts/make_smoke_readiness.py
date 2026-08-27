#!/usr/bin/env python3
"""Create a synthetic readiness ledger solely for selection-path testing."""
import argparse
import hashlib
import os

import pandas as pd


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis_panel.tsv")
    parser.add_argument("--h2", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    smoke_root = os.path.normpath("results/_smoketest")
    out_path = os.path.normpath(args.out)
    if out_path != smoke_root and not out_path.startswith(smoke_root + os.sep):
        raise SystemExit("ERROR: synthetic readiness must stay under results/_smoketest/")

    panel = pd.read_csv(args.config, sep="\t", dtype=str).fillna("")
    h2 = pd.read_csv(args.h2, sep="\t", dtype=str).fillna("")
    verdict = dict(zip(h2["trait"], h2["verdict"]))
    reason = dict(zip(h2["trait"], h2["qc_reason"]))
    panel_hash = file_sha256(args.config)
    rows = []
    for _, trait in panel.iterrows():
        passed = verdict.get(trait["trait_id"]) == "PASS"
        rows.append({
            "panel_version": trait["panel_version"],
            "manifest_sha256": panel_hash,
            "trait_id": trait["trait_id"],
            "label": trait["label"],
            "domain": trait["domain"],
            "type": trait["type"],
            "declared_source_status": "SYNTHETIC",
            "source_verified": True,
            "harmonization_ready": True,
            "ldsc_ready": True,
            "liability_h2_ready": True,
            "phase1_pass": passed,
            "readiness_stage": "PHASE1_PASS" if passed else "LIABILITY_H2_READY",
            "h2_verdict": verdict.get(trait["trait_id"], "UNASSESSED"),
            "h2_qc_reason": reason.get(trait["trait_id"], "UNASSESSED"),
        })
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Wrote SYNTHETIC readiness ledger: {args.out}")


if __name__ == "__main__":
    main()
