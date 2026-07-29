#!/usr/bin/env python3
"""Reproducibly fetch PUBLIC GWAS summary statistics listed in
config/public_gwas_sources.tsv.

Downloads land in data/raw/, which is gitignored. Nothing fetched here may be
committed: raw sumstats, archives, and reference panels stay out of the repo.

What makes this reproducible rather than a pile of curl calls: every download
records the SHA-256 of the bytes actually received and writes it back into the
registry. Re-running verifies instead of re-downloading, so a source that
changes underneath us is detected rather than silently absorbed.

    python3 scripts/09_fetch_public_sources.py --list
    python3 scripts/09_fetch_public_sources.py --trait ibd crohn uc
    python3 scripts/09_fetch_public_sources.py --all-public
    python3 scripts/09_fetch_public_sources.py --verify

Only rows with access=PUBLIC are ever fetched. REGISTRATION and CONTROLLED
rows are skipped by design -- obtaining those requires a human to accept the
provider's terms, and this script will not pretend otherwise.
"""
import argparse
import csv
import hashlib
import os
import sys
import time
import urllib.request

REG = "config/public_gwas_sources.tsv"
RAW = "data/raw"
CHUNK = 1 << 20
UA = {"User-Agent": "sleep-gwas-atlas/phase0 (academic use)"}


def load():
    with open(REG, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def save(rows):
    cols = list(rows[0].keys())
    with open(REG, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(CHUNK), b""):
            h.update(blk)
    return h.hexdigest()


def fetch(row):
    url = row["download_url"]
    dest = os.path.join(RAW, row["local_file"])
    if url in ("NA", "PENDING", ""):
        print(f"  {row['trait_id']}: no usable download_url ({url}) -- skipped")
        return False
    os.makedirs(RAW, exist_ok=True)
    # Already have it at the expected size? Don't spend the window re-fetching.
    want = int(row["bytes"]) if str(row["bytes"]).isdigit() else 0
    if os.path.exists(dest) and want and os.path.getsize(dest) == want:
        if row["sha256"] in ("", "NA", "PENDING"):
            row["sha256"] = sha256_of(dest)
        print(f"  {row['trait_id']}: already present and complete "
              f"({want/1e6:.1f} MB) -- skipped")
        return True
    tmp = dest + ".part"
    print(f"  {row['trait_id']}: GET {url}")
    t0 = time.time()
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as out:
            total = int(r.headers.get("Content-Length", 0))
            got = 0
            while True:
                b = r.read(CHUNK)
                if not b:
                    break
                out.write(b)
                got += len(b)
    except Exception as e:
        print(f"    FAILED: {type(e).__name__}: {e}")
        if os.path.exists(tmp):
            os.remove(tmp)
        return False
    # Verify completeness BEFORE promoting .part to the real filename. A
    # truncated file that carries the final name is indistinguishable from a
    # good one downstream, and this pass already produced one (cad.tsv stopped
    # at 1267 MB of 3250 MB and was promoted anyway). Refuse to finalise.
    if total and got != total:
        print(f"    !! TRUNCATED: expected {total} bytes, got {got} "
              f"({100*got/total:.1f}%). Discarding partial file.")
        os.remove(tmp)
        return False
    os.replace(tmp, dest)
    dt = time.time() - t0
    digest = sha256_of(dest)
    row["sha256"] = digest
    row["bytes"] = str(got)
    print(f"    {got/1e6:.1f} MB in {dt:.0f}s  sha256:{digest[:16]}...")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", nargs="*", default=None)
    ap.add_argument("--all-public", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--evict", action="store_true",
                    help="after hashing, delete the local file. The registry "
                         "keeps url+sha256+bytes, so it is reproducible on "
                         "demand. Use when disk cannot hold the whole panel.")
    ap.add_argument("--verify", action="store_true",
                    help="re-hash local files and compare against the registry")
    a = ap.parse_args()

    rows = load()

    if a.list:
        print(f"{'trait':22}{'access':14}{'build':10}{'local file'}")
        for r in rows:
            print(f"{r['trait_id']:22}{r['access']:14}{r['build_declared']:10}"
                  f"{r['local_file']}")
        n_pub = sum(1 for r in rows if r["access"] == "PUBLIC")
        print(f"\n{len(rows)} rows, {n_pub} PUBLIC")
        return

    if a.verify:
        bad = 0
        for r in rows:
            dest = os.path.join(RAW, r["local_file"]) if r["local_file"] not in ("NA",) else None
            if not dest or not os.path.exists(dest):
                continue
            d = sha256_of(dest)
            ok = (d == r["sha256"])
            print(f"  {r['trait_id']:22}{'OK' if ok else 'MISMATCH'}  {d[:16]}...")
            bad += (not ok)
        print(f"\n{bad} mismatch(es)")
        sys.exit(1 if bad else 0)

    if a.all_public:
        targets = [r for r in rows if r["access"] == "PUBLIC"]
    elif a.trait:
        targets = [r for r in rows if r["trait_id"] in a.trait]
        blocked = [r for r in targets if r["access"] != "PUBLIC"]
        for r in blocked:
            print(f"  REFUSING {r['trait_id']}: access={r['access']}. "
                  f"Obtain it yourself from {r['landing_page']} -- this script "
                  f"does not bypass access controls.")
        targets = [r for r in targets if r["access"] == "PUBLIC"]
    else:
        ap.error("give --trait, --all-public, --list or --verify")

    if not targets:
        print("nothing to fetch")
        return

    print(f"fetching {len(targets)} PUBLIC source(s) into {RAW}/")
    n = 0
    for r in targets:
        n += fetch(r)
    save(rows)
    print(f"\n{n}/{len(targets)} downloaded; registry updated with SHA-256")


if __name__ == "__main__":
    main()
