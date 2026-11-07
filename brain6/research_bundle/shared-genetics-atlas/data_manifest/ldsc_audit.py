#!/usr/bin/env python3
"""Check source-byte receipts for the two published global re-executions."""
import csv, hashlib, json
from pathlib import Path
root=Path(__file__).resolve().parent.parent
out=Path(__file__).resolve().parent
receipts={name:root/'global_audit'/name/'receipt.json' for name in ('insomnia_adhd_rerun','fi_ffs_rerun')}
expected={};links={}
for name,receipt in receipts.items():
    for path,h in json.loads(receipt.read_text())['source_sha256'].items():
        p=Path(path)
        include=(path.startswith('/Volumes/Extreme SSD/') and
                 (p.name in {'insomnia.sumstats.gz','adhd.sumstats.gz','frailty.sumstats.gz','fried_frailty_score.sumstats.gz','ldsc.py'}
                  or p.parent.name=='eur_w_ld_chr' and (p.name.endswith('.l2.ldscore.gz') or p.name.endswith('.l2.M_5_50'))
                  or p.parent.name=='ldscore' and p.suffix=='.py' and not p.name.startswith('._')))
        if not include:continue
        if path in expected:assert expected[path]==h
        expected[path]=h;links.setdefault(path,[]).append(name)
assert len(expected)>=50
rows=[]
for path in sorted(expected):
    p=Path(path);status='missing';digest='';size=''
    if p.is_file():
        h=hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
        digest=h.hexdigest();size=p.stat().st_size;status='match' if digest==expected[path] else 'changed'
    role=('munged GWAS' if p.name.endswith('.sumstats.gz') else 'LD reference' if p.parent.name=='eur_w_ld_chr' else 'LDSC source code')
    rows.append({'role':role,'path':path,'bytes':size,'expected_sha256':expected[path],
                 'observed_sha256':digest,'status':status,'source_receipts':','.join(links[path])})
with (out/'ldsc_inputs.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
summary={'file_count':len(rows),'counts':{s:sum(r['status']==s for r in rows) for s in ('match','changed','missing')},
 'role_counts':{s:sum(r['role']==s for r in rows) for s in ('munged GWAS','LD reference','LDSC source code')},
 'scope':'Byte check of receipt-selected inputs only; not upstream origin, a software environment lock, scientific validity or clinical use.'}
(out/'ldsc_check.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if summary['counts']['changed'] or summary['counts']['missing']:raise SystemExit(1)
