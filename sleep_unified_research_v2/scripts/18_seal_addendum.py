#!/usr/bin/env python3
"""Seal only small research evidence files, never SSD GWAS bodies."""
import argparse
import hashlib
import json
from pathlib import Path

PACKAGE=Path(__file__).resolve().parents[1]
MANIFEST=PACKAGE/'hashes/ARTIFACT_SHA256_v2.tsv'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def files():
    return sorted(p for p in PACKAGE.rglob('*') if p.is_file() and p!=MANIFEST
                  and '__pycache__' not in p.parts and not p.name.startswith('._'))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.verify:
        rows=[line.split('\t') for line in MANIFEST.read_text().splitlines()[1:]]
        for expected,n,relative in rows:
            path=PACKAGE/relative
            if path.stat().st_size!=int(n) or sha(path)!=expected:
                raise SystemExit('ARTIFACT_HASH_FAILED:'+relative)
        listed={r[2] for r in rows}
        actual={str(p.relative_to(PACKAGE)) for p in files()}
        if listed!=actual:
            raise SystemExit('MANIFEST_CARDINALITY_FAILED')
        print('ADDENDUM_SHA256_PASS files='+str(len(rows)))
        return
    if MANIFEST.exists():
        raise SystemExit('PRIOR_SEAL_PRESERVED_NO_OVERWRITE')
    paths=files()
    if any(p.stat().st_size>5*1024**2 for p in paths):
        raise SystemExit('UNEXPECTED_LARGE_ARTIFACT_REQUIRES_REVIEW')
    MANIFEST.parent.mkdir(parents=True,exist_ok=True)
    with MANIFEST.open('x') as f:
        f.write('sha256\tbytes\tpath\n')
        for p in paths:
            f.write(f'{sha(p)}\t{p.stat().st_size}\t{p.relative_to(PACKAGE)}\n')
    print(json.dumps({'sealed_files':len(paths),'manifest':str(MANIFEST),'manifest_sha256':sha(MANIFEST)}))

if __name__=='__main__':
    main()
