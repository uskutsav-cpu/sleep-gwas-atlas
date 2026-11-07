#!/usr/bin/env python3
import collections,csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=p.parents[4]
variants=atlas/'variant_context/rs62264764_dense_qtl/multisignal_pilot/variants.tsv'; brain=atlas/'variant_context/rs62264764_brain_xqtl/ldz_readiness/ldz_details.tsv'; gwas=atlas/'variant_context/rs62264764_dense_qtl/multisignal_pilot/kriging_diagnostics.tsv'
def read(q):
 with q.open() as f:return list(csv.DictReader(f,delimiter='\t'))
v=read(variants);by={x['variant']:x for x in v};reasons=collections.defaultdict(list)
for path,column,contexts in [(brain,'context',['oli','exc']),(gwas,'fit',['insomnia','adhd'])]:
 rows=read(path)
 for ctx in contexts:
  subset=sorted((x for x in rows if x[column]==ctx),key=lambda x:(-abs(float(x['z_std_diff'])),x['variant']))[:10]
  for x in subset:reasons[x['variant']].append(ctx+'_top10_abs_residual')
for pos in ('117920728','117923158'):
 matches=[x for x in v if x['pos38']==pos];assert len(matches)==1;reasons[matches[0]['variant']].append('lead')
selected=[]
for name in sorted(reasons,key=lambda z:int(by[z]['pos38'])):
 x=by[name];selected.append({'variant':name,'pos38':x['pos38'],'ref':x['ref'],'alt':x['alt'],'reasons':';'.join(reasons[name]),'reference_index_0':v.index(x)})
assert len({x['pos38'] for x in selected})==len(selected)
with (p/'selected.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=selected[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
(p/'patterns.txt').write_text(''.join('\t'+x['pos38']+'\t\n' for x in selected))
(p/'selection_receipt.json').write_text(json.dumps({'n':len(selected),'input_sha256':{str(q.relative_to(atlas)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (variants,brain,gwas)},'rule':'union of top ten absolute residuals in four fits plus two leads','selection_is_posthoc':True},indent=2)+'\n');print('selected',len(selected),'variants')
