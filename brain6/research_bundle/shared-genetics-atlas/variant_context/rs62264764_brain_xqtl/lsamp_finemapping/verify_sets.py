import csv,collections,json,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'study_credible_set_rows.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
assert len(r)==13 and len({(x['cell'],x['variant'],x['cs95']) for x in r})==13
assert all(x['target_id']=='ENSG00000185565' and 0<float(x['pip'])<=1 and 'is_hmt_signif=true' in x['qtl_metadata'] for x in r)
groups=collections.defaultdict(list)
for x in r:groups[x['cell']].append(x)
assert {k:len(v) for k,v in groups.items()}=={'oli':10,'exc':3}
for x in groups.values():assert len({r['cs95'] for r in x})==1 and sum(float(r['pip']) for r in x)>=.95
assert all(x['variant'] not in ('chr3:117920728:G:A','chr3:117923158:A:C') for x in r)
old=p.parents[1]/'rs62264764_dense_qtl/multisignal_pilot/credible_sets.tsv'
with old.open() as f:prior=list(csv.DictReader(f,delimiter='\t'))
sets={x['fit']:{s.replace('_',':') for s in x['variants'].split(',')} for x in prior};overlap={}
for cell,d in groups.items():
 x={z['variant'] for z in d}
 for fit,y in sets.items():overlap[cell+'_'+fit]=sorted(x&y)
out={'rows':len(r),'counts':{c:len(d) for c,d in groups.items()},'pip_sums':{c:sum(float(x['pip']) for x in d) for c,d in groups.items()},'rs62264764_and_rs62264767_in_selected_sets':False,'overlaps_with_unvalidated_prior_susie_sets':overlap,'warning':'PIP sums and simple set intersection are not colocalization or proof of variant causality; prior disease fits have LD–Z concerns.'}
(p/'set_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
