#!/usr/bin/env python3
"""Globally count duplicate coordinate/allele keys in latent Catalog GWAS files.

Uses external sort temporary storage on the configured data volume. Source files
are read-only; no variants are dropped or rewritten.
"""
from __future__ import annotations
import argparse, csv, os, shutil, subprocess, tempfile
from pathlib import Path

FIELDS=['accession','file','rows_sorted','duplicate_key_groups','duplicate_excess_rows','global_key_uniqueness','key_definition','method']
AWK='BEGIN { FS="\\t"; OFS="\\t" } NR>1 { print $1, $2, toupper($3), toupper($4) }'

def audit(path:Path, accession:str, external_root:Path):
    temp=Path(tempfile.mkdtemp(prefix='latent-variant-sort-',dir=external_root))
    awk=subprocess.Popen(['awk',AWK,str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    sort=subprocess.Popen(['sort','-S','1G','-T',str(temp),'-t','\t','-k1,1','-k2,2n','-k3,3','-k4,4'],stdin=awk.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env={**os.environ,'LC_ALL':'C'})
    assert awk.stdout is not None and sort.stdout is not None
    awk.stdout.close()
    rows=groups=excess=0; previous=None
    for raw in sort.stdout:
        key=raw.rstrip(b'\n'); rows+=1
        if key==previous:
            excess+=1
            if excess==1 or previous_group!=key: groups+=1
        previous_group=key
        previous=key
    sort.stdout.close()
    sort_error=sort.stderr.read() if sort.stderr else b''
    awk_error=awk.stderr.read() if awk.stderr else b''
    sort_code=sort.wait(); awk_code=awk.wait()
    shutil.rmtree(temp)
    if sort_code or awk_code:
        raise RuntimeError(f'{accession}: awk={awk_code} sort={sort_code}; {awk_error.decode(errors="replace")} {sort_error.decode(errors="replace")}')
    return {'accession':accession,'file':str(path),'rows_sorted':rows,'duplicate_key_groups':groups,'duplicate_excess_rows':excess,'global_key_uniqueness':'UNIQUE' if excess==0 else 'DUPLICATE_KEYS_PRESENT','key_definition':'chromosome + numeric position + uppercase effect allele + uppercase other allele','method':'awk fields 1-4 -> LC_ALL=C external sort -> adjacent duplicate count; no filtering'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]); ap.add_argument('--external-storage-root',type=Path,required=True); ap.add_argument('--output',type=Path,default=Path('frailty_paper/manifests/gwas_latent_duplicate_audit.tsv')); args=ap.parse_args(); repo=args.repo.resolve(); ext=args.external_storage_root.resolve()
    if shutil.disk_usage(ext).free < 20*1024**3: raise SystemExit('external storage has less than 20 GiB free; no audit started')
    index=list(csv.DictReader((repo/'frailty_paper/manifests/latent_frailty_accession_index.tsv').open(),delimiter='\t')); rows=[]
    for m in index:
        p=repo/'frailty_paper/data/gwas/latent_frailty_catalog'/m['accession']/f"{m['accession']}.tsv"
        result=audit(p,m['accession'],ext); result['file']=p.relative_to(repo).as_posix(); rows.append(result)
        print(f"{result['accession']} rows={result['rows_sorted']} duplicate_key_groups={result['duplicate_key_groups']} duplicate_excess_rows={result['duplicate_excess_rows']}",flush=True)
    dest=args.output if args.output.is_absolute() else repo/args.output; dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    print(f'LATENT_DUPLICATE_AUDIT_WRITTEN rows={len(rows)} output={dest}')
if __name__=='__main__':main()
