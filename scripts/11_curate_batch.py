#!/usr/bin/env python3
"""Download -> verify -> record -> (optionally) evict, one trait at a time.

Phase 0 asks that a source was downloaded, its header inspected, its build
proven from the file, and its sample handling documented. All of that leaves
evidence: a SHA-256, a byte count, a build verdict and a variant count. The
raw file itself does not have to survive -- with the URL and hash recorded the
exact bytes are reproducible on demand, and the hash proves you got them back.

That matters because disk is finite and a 100-trait panel is not going to fit.
--evict trades local storage for reproducibility, and refuses to evict
anything it could not first verify.

    python3 scripts/11_curate_batch.py --trait celiac migraine --evict
    python3 scripts/11_curate_batch.py --pending --evict --budget-gb 4

A trait is only marked verified if the build resolves to GRCh37 AND the file
has at least MIN_VARIANTS rows. Anything else is left TODO with the reason
written into the registry.
"""
import argparse
import csv
import gzip
import importlib.util
import os
import sys

REG = "config/public_gwas_sources.tsv"
TRAITS = "config/traits.tsv"
RAW = "data/raw"
MIN_VARIANTS = 1_000_000

# Importing these modules must not disturb our own command line: the audit
# module inspects sys.argv at import time, so save and restore it.
_argv = list(sys.argv)
spec = importlib.util.spec_from_file_location("audit", "scripts/10_phase0_audit.py")
audit = importlib.util.module_from_spec(spec)
sys.argv = ["audit"]
spec.loader.exec_module(audit)

fspec = importlib.util.spec_from_file_location("fetch", "scripts/09_fetch_public_sources.py")
fetcher = importlib.util.module_from_spec(fspec)
fspec.loader.exec_module(fetcher)
sys.argv = _argv


def read(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def write(path, rows):
    with open(path, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)


def count_variants(path):
    op = gzip.open if path.endswith(".gz") else open
    try:
        with op(path, "rt", errors="replace") as fh:
            return sum(1 for _ in fh)
    except Exception:
        return -1


def _promote(done):
    """Mark the given trait_ids CURATED in traits.tsv, idempotently."""
    if not done:
        return
    trows = read(TRAITS)
    j = {c: n for n, c in enumerate(trows[0])}
    for r in trows[1:]:
        if r[j["trait_id"]] in done:
            r[j["status"]] = "CURATED"
            r[j["build"]] = "hg19"
            if r[j["type"]] == "binary" and r[j["pop_prev"]] not in ("NA", "UNKNOWN"):
                r[j["pop_prev"]] = "UNKNOWN"
    write(TRAITS, trows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", nargs="*", default=None)
    ap.add_argument("--pending", action="store_true",
                    help="every PUBLIC row whose build_verified is PENDING")
    ap.add_argument("--evict", action="store_true")
    ap.add_argument("--reconcile", action="store_true",
                    help="promote from registry evidence alone. A killed batch "
                         "still wrote build_verified and n_variants per row, so "
                         "the promotions can always be rebuilt without "
                         "re-downloading anything.")
    ap.add_argument("--budget-gb", type=float, default=None,
                    help="stop once this much has been downloaded this run")
    a = ap.parse_args()

    rows = read(REG)
    h = rows[0]
    i = {c: n for n, c in enumerate(h)}

    if a.reconcile:
        good = [r[i["trait_id"]] for r in rows[1:]
                if r[i["build_verified"]] == "GRCh37"
                and str(r[i["n_variants"]]).isdigit()
                and int(r[i["n_variants"]]) >= MIN_VARIANTS]
        _promote(good)
        print(f"reconciled: {len(good)} traits promoted from registry evidence")
        return

    if a.pending:
        targets = [r for r in rows[1:]
                   if r[i["access"]] == "PUBLIC"
                   and r[i["build_verified"]] in ("PENDING", "")]
    else:
        targets = [r for r in rows[1:] if r[i["trait_id"]] in (a.trait or [])]

    if not targets:
        print("nothing to do")
        return

    done, spent = [], 0.0
    for r in targets:
        tid = r[i["trait_id"]]
        if a.budget_gb and spent >= a.budget_gb:
            print(f"  budget of {a.budget_gb} GB reached -- stopping")
            break
        dest = os.path.join(RAW, r[i["local_file"]])
        d = {k: r[i[k]] for k in i}          # fetch() takes a dict-like row
        if not os.path.exists(dest):
            if not fetcher.fetch(d):
                r[i["notes"]] = "DOWNLOAD FAILED. " + r[i["notes"]]
                continue
            r[i["sha256"]] = d["sha256"]
            r[i["bytes"]] = d["bytes"]
        spent += int(r[i["bytes"]] or 0) / 1e9

        (verdict, _, _), member = audit.check_build(dest)
        nv = count_variants(dest)
        r[i["build_verified"]] = "GRCh37" if verdict.startswith("GRCh37") else (
            "GRCh38" if verdict.startswith("GRCh38") else "UNVERIFIED")
        if "n_variants" in i:
            r[i["n_variants"]] = str(nv)
        if member and "archive_member" in i:
            r[i["archive_member"]] = member

        ok = verdict.startswith("GRCh37") and nv >= MIN_VARIANTS
        tag = "OK " if ok else "NO "
        print(f"  {tag}{tid:24}{verdict[:46]:46}{nv:>12,} variants")
        if not ok:
            reason = ("build is not GRCh37" if not verdict.startswith("GRCh37")
                      else f"only {nv:,} variants -- not genome-wide")
            r[i["notes"]] = f"NOT CURATED: {reason}. " + r[i["notes"]]
            r[i["access"]] = "UNAVAILABLE" if not verdict.startswith("GRCh37") else r[i["access"]]
        else:
            done.append(tid)
            if a.evict:
                os.remove(dest)
                r[i["notes"]] = ("EVICTED after verification: file deleted to free "
                                 "disk. url+sha256+bytes recorded, so it is exactly "
                                 "reproducible. " + r[i["notes"]])
        write(REG, rows)
        # Promote incrementally. Writing traits.tsv only at the end means a
        # killed run loses every promotion it had already earned -- and long
        # batches do get killed.
        _promote(done)

    if done:
        trows = read(TRAITS)
        j = {c: n for n, c in enumerate(trows[0])}
        for r in trows[1:]:
            if r[j["trait_id"]] in done:
                r[j["status"]] = "CURATED"
                r[j["build"]] = "hg19"
                if r[j["type"]] == "binary" and r[j["pop_prev"]] not in ("NA", "UNKNOWN"):
                    r[j["pop_prev"]] = "UNKNOWN"
        write(TRAITS, trows)
    print(f"\ncurated {len(done)}: {' '.join(done)}")


if __name__ == "__main__":
    main()
