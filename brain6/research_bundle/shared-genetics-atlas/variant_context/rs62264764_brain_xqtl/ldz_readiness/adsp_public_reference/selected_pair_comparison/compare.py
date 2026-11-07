#!/usr/bin/env python3
import collections,csv,gzip,hashlib,itertools,json,math,statistics
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;atlas=p.parents[4];prior=atlas/'variant_context/rs62264764_dense_qtl/multisignal_pilot';m=prior/'ld_matrix.tsv.gz';v=prior/'variants.tsv';receipt=json.loads((prior/'ld_receipt.json').read_text())
def digest(q):return hashlib.sha256(q.read_bytes()).hexdigest()
def read(q):
 with q.open() as f:return list(csv.DictReader(f,delimiter='\t'))
assert digest(m)==receipt['ld_sha256'] and digest(v)==receipt['variants_sha256']
R=np.loadtxt(m);assert R.shape==(1548,1548) and np.max(np.abs(R-R.T))<1e-10
sel=read(p/'selected.tsv');variants=read(v)
for x in sel:
 a=variants[int(x['reference_index_0'])];assert a['variant']==x['variant'] and a['pos38']==x['pos38'] and a['ref']==x['ref'] and a['alt']==x['alt']
raw=read(p/'source_pairs.tsv');groups=collections.defaultdict(list)
for x in raw:groups[tuple(sorted((x['POS1'],x['POS2']),key=int))].append(x)
out=[];c=collections.Counter()
for a,b in itertools.combinations(sel,2):
 key=tuple(sorted((a['pos38'],b['pos38']),key=int));matches=[]
 for x in groups.get(key,[]):
  alleles={x['POS1']:x['REF:ALT1'],x['POS2']:x['REF:ALT2']}
  if alleles.get(a['pos38'])==a['ref']+':'+a['alt'] and alleles.get(b['pos38'])==b['ref']+':'+b['alt']:matches.append(x)
 status='matched' if len(matches)==1 else 'no_exact_record' if not matches else 'duplicate_exact_record';c[status]+=1
 x={'variant1':a['variant'],'variant2':b['variant'],'status':status,'kgp_r2':float(R[int(a['reference_index_0']),int(b['reference_index_0'])]**2),'adsp_r2':'','absolute_difference':''}
 if status=='matched':
  f=float(matches[0]['Rsq']);assert math.isfinite(f) and 0<=f<=1;x['adsp_r2']=f;x['absolute_difference']=abs(f-x['kgp_r2'])
 out.append(x)
with (p/'pair_comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
d=sorted(float(x['absolute_difference']) for x in out if x['status']=='matched');assert d
s={'variants':len(sel),'all_possible_pairs':len(out),'status_counts':dict(c),'matched_abs_r2_difference_median':statistics.median(d),'matched_abs_r2_difference_p95':d[int(.95*(len(d)-1))],'matched_abs_r2_difference_max':max(d),'matched_difference_ge_0.1':sum(x>=.1 for x in d),'matched_difference_ge_0.2':sum(x>=.2 for x in d),'source_rows':len(raw),'input_sha256':{x.name:digest(x) for x in (m,v,p/'selected.tsv',p/'source_pairs.tsv')},'interpretation':'Selected returned pairs only; missing is not zero; sample and LD-statistic definitions differ.'}
(p/'comparison_receipt.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
