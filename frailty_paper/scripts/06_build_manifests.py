#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, os
from pathlib import Path

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo', required=True); args=ap.parse_args()
    root=Path(args.repo).resolve(); paper=root/'frailty_paper'
    rows=[]
    for logical_base in [root/'data/raw', paper/'data']:
        if not logical_base.exists(): continue
        physical_base=logical_base.resolve()
        if not physical_base.is_dir(): continue
        for p in sorted(physical_base.rglob('*')):
            if p.is_file() and p.name != '.DS_Store' and not p.name.startswith('._'):
                rel=logical_base.relative_to(root)/p.relative_to(physical_base)
                rows.append({'path':str(rel),'bytes':p.stat().st_size,'sha256':sha256(p)})
    out=paper/'manifests/file_checksums.tsv'; out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t'); w.writeheader(); w.writerows(rows)

    targets={'frailty','grip_strength','parental_lifespan','longevity','cognitive_performance','alz','parkinson',
             'insomnia','sleepdur','shortsleep','longsleep','chronotype','sleepiness','napping','snoring','sleep_apnea','sleep_efficiency','accel_sleep_duration','sleep_timing'}
    trait_rows=[]
    traits=root/'config/traits.tsv'
    if traits.exists():
        with traits.open(encoding='utf-8',newline='') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r.get('trait_id') in targets:
                    raw=root/'data/raw'/r.get('raw_file','')
                    trait_rows.append({
                        'trait_id':r.get('trait_id',''),'label':r.get('label',''),'domain':r.get('domain',''),
                        'source_note':r.get('source_note',''),'raw_file':r.get('raw_file',''),'n_total':r.get('n_total',''),
                        'build':r.get('build',''),'pmid':r.get('pmid',''),'ancestry':r.get('ancestry',''),
                        'raw_present':'YES' if raw.is_file() and raw.stat().st_size>0 else 'NO',
                        'raw_bytes':raw.stat().st_size if raw.is_file() else 0,
                    })
    out2=paper/'manifests/existing_repo_trait_audit.tsv'
    fields=['trait_id','label','domain','source_note','raw_file','n_total','build','pmid','ancestry','raw_present','raw_bytes']
    with out2.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(trait_rows)
    print('Wrote',out); print('Wrote',out2)
if __name__=='__main__': main()
