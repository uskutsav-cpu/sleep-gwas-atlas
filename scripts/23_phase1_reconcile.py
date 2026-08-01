#!/usr/bin/env python3
"""Reconcile Phase 1 counts across filesystem, ledger and summary tables.

This exists because a real inconsistency shipped: PROJECT_STATE reported 15
traits complete while h2_summary.tsv held 14 rows and recorded 0 failures. The
cause was not a miscount but STALENESS -- summary tables are snapshots, and the
batch kept completing traits after they were written. Nothing detected the
drift because nothing compared them.

Exit status is non-zero on any mismatch, so it can gate a milestone commit.

    python3 scripts/23_phase1_reconcile.py            # check
    python3 scripts/23_phase1_reconcile.py --repair   # regenerate then check
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys

LEDGER = "results/phase_ledger.tsv"
H2 = "results/phase1/h2_summary.tsv"
H2_PAT = re.compile(r"Total (Observed|Liability) scale h2:\s*(-?[\d.eE+-]+)")


def fs_state():
    munged = {os.path.basename(p).replace(".sumstats.gz", "")
              for p in glob.glob("data/munged/*.sumstats.gz")}
    h2_ok, h2_bad = set(), set()
    for p in glob.glob("results/logs/h2_*.log"):
        t = os.path.basename(p)[3:-4]
        if H2_PAT.search(open(p, errors="replace").read()):
            h2_ok.add(t)
        else:
            h2_bad.add(t)
    return munged, h2_ok, h2_bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repair", action="store_true")
    a = ap.parse_args()

    if a.repair:
        subprocess.run([".ldsc-env/bin/python", "scripts/05_collate.py",
                        "--mode", "h2", "--logdir", "results/logs",
                        "--out", H2], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        subprocess.run([".ldsc-env/bin/python", "scripts/13_build_ledger.py"],
                       stdout=subprocess.DEVNULL)

    munged, h2_ok, h2_bad = fs_state()
    complete = munged & h2_ok
    problems = []

    if not os.path.exists(H2):
        problems.append(f"{H2} missing")
        summ = set()
    else:
        rows = list(csv.DictReader(open(H2), delimiter="\t"))
        summ = {r["trait"] for r in rows}
        npass = sum(1 for r in rows if r.get("verdict") == "PASS")
        nfail = sum(1 for r in rows if r.get("verdict") == "DROP")
        if len(rows) != npass + nfail:
            problems.append(f"h2_summary: {len(rows)} rows but "
                            f"{npass} PASS + {nfail} DROP")
        print(f"h2_summary rows={len(rows)}  PASS={npass}  DROP={nfail}")

    print(f"filesystem: munged={len(munged)} h2_parsed={len(h2_ok)} "
          f"complete={len(complete)}")

    if complete - summ:
        problems.append(f"complete on disk but ABSENT from h2_summary: "
                        f"{sorted(complete - summ)}")
    if summ - complete:
        problems.append(f"in h2_summary but NOT complete on disk: "
                        f"{sorted(summ - complete)}")
    if h2_bad:
        problems.append(f"h2 logs with no parsable h2 (likely killed "
                        f"mid-write): {sorted(h2_bad)}")

    if os.path.exists(LEDGER):
        led = list(csv.DictReader(open(LEDGER), delimiter="\t"))
        led_done = {r["trait_id"] for r in led if r["h2_status"] == "DONE"}
        print(f"ledger: {len(led)} traits, h2_status=DONE for {len(led_done)}")
        if led_done != complete:
            problems.append(f"ledger DONE != filesystem complete; "
                            f"ledger-only={sorted(led_done - complete)} "
                            f"disk-only={sorted(complete - led_done)}")

    print()
    if problems:
        print(f"RECONCILIATION FAILED ({len(problems)} problem(s)):")
        for p in problems:
            print(f"  - {p}")
        sys.exit(1)
    print(f"RECONCILED: {len(complete)} traits complete, counts agree across "
          f"filesystem, ledger and h2_summary")


if __name__ == "__main__":
    main()
