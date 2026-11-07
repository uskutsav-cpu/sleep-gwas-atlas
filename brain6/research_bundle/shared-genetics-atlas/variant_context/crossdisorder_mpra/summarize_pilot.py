import csv,json,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent;ssd=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra')
with (p/'ranked_coverage.tsv').open() as f: coverage=list(csv.DictReader(f,delimiter='\t'))
models={}
for m in ['v124','v134']:
 with (ssd/f'limma_pilot_{m}_all.tsv').open() as f: models[m]={r['variant']:r for r in csv.DictReader(f,delimiter='\t')}
rows=[]
for x in coverage:
 key=x['source_key']
 if not all(key in models[m] for m in models):continue
 r={k:x[k] for k in ('rsid','region','grch38','source_key','listed_excluded')}
 for m in models:
  for c in ['logFC','P.Value','adj.P.Val']:r[m+'_'+c]=models[m][key][c]
 rows.append(r)
assert len(rows)==42
with (p/'provisional_allelic_42.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
counts={m:sum(float(r[m+'_adj.P.Val'])<.05 for r in rows) for m in models}; flips=sum((float(r['v124_adj.P.Val'])<.05)!=(float(r['v134_adj.P.Val'])<.05) for r in rows)
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
receipt={'count':len(rows),'fdr_under_0.05_by_mode':counts,'changed_threshold_status':flips,'models':{m:{'all_rows':len(models[m]),'all_fdr_under_0.05':sum(float(v['adj.P.Val'])<.05 for v in models[m].values()),'output_sha256':sha(ssd/f'limma_pilot_{m}_all.tsv')} for m in models},'coverage_sha256':sha(p/'ranked_coverage.tsv'),'pilot_code_sha256':sha(p/'limma_pilot.R'),'output_sha256':sha(p/'provisional_allelic_42.tsv')}
(p/'pilot_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print('rs4916723',next(r for r in rows if r['rsid']=='rs4916723'))
