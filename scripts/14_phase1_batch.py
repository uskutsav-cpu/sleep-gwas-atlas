#!/usr/bin/env python3
"""Resumable per-trait Phase 1 pipeline: fetch -> harmonize -> munge -> h2 -> evict.

Finishing Phase 1 for ~90 traits is not one long command. Each trait costs a
few hundred MB of download, several minutes of chunked harmonization, and this
machine has 8 GB of RAM and single-digit GB of free disk. So the unit of work
is ONE TRAIT, fully completed and checkpointed, and the driver is safe to kill
and restart at any point.

Resume logic is state-based, not position-based: a trait is skipped when its
munged sumstats and h2 log already exist. Nothing is recomputed to "catch up".

Disk discipline: raw and harmonized files are deleted once the munged file
exists, because munged sumstats are ~8 MB while raw+harmonized are ~700 MB.
The registry keeps url + sha256 + bytes, so any input is exactly reproducible
(verified in practice -- re-downloads have returned byte-identical hashes).

    python3 scripts/14_phase1_batch.py --status
    python3 scripts/14_phase1_batch.py --role core_sleep disease --limit 6
    python3 scripts/14_phase1_batch.py --trait ibd ra --keep-raw

Every trait ends in an explicit terminal state written to the ledger; a
failure never silently vanishes.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import time

REG = "config/public_gwas_sources.tsv"
TRAITS = "config/traits.tsv"
LEDGER = "results/phase_ledger.tsv"
PY = ".ldsc-env/bin/python"
RAW, HARM, MUNGED, LOGS = "data/raw", "data/harmonized", "data/munged", "results/logs"
CHUNK = 1_500_000


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(cmd, log, timeout=3000):
    """Run a step, capturing everything. Never raises; returns (ok, note)."""
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "w") as fh:
        try:
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                               timeout=timeout)
        except subprocess.TimeoutExpired:
            return False, f"TIMEOUT after {timeout}s"
    return (p.returncode == 0), f"exit={p.returncode}"


def done(tid):
    return (os.path.exists(f"{MUNGED}/{tid}.sumstats.gz")
            and os.path.exists(f"{LOGS}/h2_{tid}.log"))


def evict(tid, src_row, keep_raw):
    freed = 0
    paths = [f"{HARM}/{tid}.harmonized.tsv.gz"]
    if not keep_raw:
        lf = src_row.get("local_file", "")
        if lf and lf != "NA":
            paths.append(os.path.join(RAW, lf))
    for p in paths:
        if os.path.exists(p):
            freed += os.path.getsize(p)
            os.remove(p)
    return freed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", nargs="*", default=None)
    ap.add_argument("--role", nargs="*", default=None,
                    help="core_sleep / secondary_sleep / disease / aging")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--min-free-gb", type=float, default=2.0,
                    help="stop before disk drops below this")
    a = ap.parse_args()

    src = {r["trait_id"]: r for r in read(REG)}
    ledger = {r["trait_id"]: r for r in read(LEDGER)} if os.path.exists(LEDGER) else {}

    if a.status:
        tot = len(ledger)
        d = sum(1 for t in ledger if done(t))
        print(f"ledger traits: {tot}   phase1 complete (munged+h2): {d}")
        for role in ("core_sleep", "secondary_sleep", "disease", "aging"):
            ids = [t for t, r in ledger.items() if r["role"] == role]
            print(f"  {role:16} {sum(1 for t in ids if done(t)):3}/{len(ids)}")
        return

    if a.trait:
        targets = list(a.trait)
    else:
        roles = set(a.role or ["core_sleep", "disease"])
        # Only PUBLIC rows with a real size are candidates. Filtering here
        # rather than inside the loop matters: an UNAVAILABLE row has bytes=0,
        # so it sorted FIRST and consumed the whole --limit on skips.
        targets = [t for t, r in ledger.items()
                   if r["role"] in roles
                   and src.get(t, {}).get("access") == "PUBLIC"
                   and str(src.get(t, {}).get("bytes") or "0").isdigit()
                   and int(src[t]["bytes"]) > 0
                   # exclude documented blockers HERE, not in the loop: they
                   # are still access=PUBLIC, so filtering later let them
                   # consume --limit on skips
                   and "PHASE1_BLOCKED" not in (src[t].get("notes") or "")]
        targets.sort(key=lambda t: int(src[t]["bytes"]))

    todo = [t for t in targets if not done(t)]
    if a.limit:
        todo = todo[:a.limit]
    if not todo:
        print("nothing to do -- every requested trait already has munged + h2")
        return
    print(f"{len(todo)} trait(s) to process")

    ok_n = fail_n = 0
    for tid in todo:
        free_gb = shutil.disk_usage(".").free / 1e9
        if free_gb < a.min_free_gb:
            print(f"  STOP: only {free_gb:.1f} GB free (< {a.min_free_gb})")
            break
        s = src.get(tid)
        if not s:
            print(f"  {tid:22} SKIP no source row"); fail_n += 1; continue
        if s["access"] != "PUBLIC":
            print(f"  {tid:22} SKIP access={s['access']}"); continue
        # A trait with a documented terminal Phase 1 blocker must not be
        # retried every pass -- the loop above burned three whole rounds
        # re-failing the same three traits.
        if "PHASE1_BLOCKED" in (s.get("notes") or ""):
            print(f"  {tid:22} SKIP phase1-blocked (see registry notes)"); continue

        t0 = time.time()
        lf = s["local_file"]
        raw = os.path.join(RAW, lf)

        if not os.path.exists(raw):
            ok, note = run([PY, "scripts/09_fetch_public_sources.py",
                            "--trait", tid], f"{LOGS}/fetch_{tid}.log", 2400)
            if not ok or not os.path.exists(raw):
                print(f"  {tid:22} FAIL fetch ({note})"); fail_n += 1; continue

        harm = f"{HARM}/{tid}.harmonized.tsv.gz"
        if not os.path.exists(harm) or os.path.getsize(harm) == 0:
            if os.path.exists(harm):
                os.remove(harm)          # 0-byte leftover from a killed run
            ok, note = run([PY, "scripts/01_harmonize.py", "--trait", tid,
                            "--config", TRAITS, "--infile", raw,
                            "--outdir", HARM, "--chunksize", str(CHUNK)],
                           f"{LOGS}/harmonize_{tid}.log", 3000)
            if not ok or not os.path.exists(harm) or os.path.getsize(harm) == 0:
                print(f"  {tid:22} FAIL harmonize ({note})"); fail_n += 1; continue

        ok, note = run(["bash", "scripts/12_phase1_run.sh", tid],
                       f"{LOGS}/phase1_{tid}.log", 3000)
        # LDSC finishes writing its .log slightly after the process returns, so
        # an immediate done() check can report a false failure (it did for dbp,
        # which had in fact succeeded). Give it a moment before deciding.
        for _ in range(10):
            if done(tid):
                break
            time.sleep(1)
        if not done(tid):
            print(f"  {tid:22} FAIL munge/h2 ({note})"); fail_n += 1; continue

        freed = evict(tid, s, a.keep_raw)
        print(f"  {tid:22} OK  {time.time()-t0:5.0f}s  freed {freed/1e6:6.0f} MB")
        ok_n += 1

    print(f"\n{ok_n} succeeded, {fail_n} failed")
    subprocess.run([PY, "scripts/13_build_ledger.py"],
                   stdout=subprocess.DEVNULL)


if __name__ == "__main__":
    main()
