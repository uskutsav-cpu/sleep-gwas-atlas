#!/usr/bin/env python3
"""Portable SHA-256 manifest; no scientific completion inference."""
import argparse
import csv
import hashlib
from pathlib import Path

PACKAGE=Path(__file__).resolve().parents[1]
MANIFEST=PACKAGE/'hashes/FILES_SHA256.tsv'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def files():
    return sorted(p for p in PACKAGE.rglob('*') if p.is_file() and p!=MANIFEST and '__pycache__' not in p.parts and not p.name.startswith('._') and p.suffix!='.pyc')
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:
        with MANIFEST.open(newline='') as f:recorded=list(csv.DictReader(f,delimiter='\t'))
        paths={r['path'] for r in recorded}
        if paths!={str(p.relative_to(PACKAGE)) for p in files()}:raise SystemExit('PACKAGE_FILESET_CHANGED')
        for r in recorded:
            p=PACKAGE/r['path']
            if not p.is_file() or p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:raise SystemExit('HASH_MISMATCH '+r['path'])
        print(f'PACKAGE_SHA256_PASS files={len(recorded)}');return
    MANIFEST.parent.mkdir(exist_ok=True)
    with MANIFEST.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader()
        for p in files():w.writerow(dict(path=str(p.relative_to(PACKAGE)),bytes=p.stat().st_size,sha256=sha(p)))
    print(f'PACKAGE_SEALED files={len(files())}')
if __name__=='__main__':main()
