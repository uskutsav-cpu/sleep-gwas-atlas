import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas'); evidence=atlas.parent/'brain6-alphagenome-pilot/results/variant_evidence.tsv'; reporter=p.parent/'comparison.tsv'
def read(path):
 with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a={r['snp']:r for r in read(evidence) if r['snp'] in ('rs325485','rs6874138')};b={r['rsid']:r for r in read(reporter) if r['rsid'] in a};assert set(a)==set(b)=={'rs325485','rs6874138'}
rows=[]
for rs in sorted(a):
 x=a[rs];y=b[rs];assert x['mapping_status']=='READY_UNIQUE_REFERENCE_VERIFIED' and x['atlas_variant']==y['variant38'];assert x['source_a1']==x['ref_grch37']==x['ref_grch38'] and x['source_a2']==x['alt_grch38'];assert float(y['published_logfc'])>0
 row={'rsid':rs,'variant_grch38':x['atlas_variant'],'conditionally_oriented_reporter_alt_minus_ref_log2':float(y['published_logfc']),'dnase_alt_minus_ref_501':float(y['dnase_501_alt_vs_ref'])}
 for trait in ('insomnia','adhd'):
  beta=float(x[trait+'_beta_a1']);se=float(x[trait+'_se']);assert math.isfinite(beta) and se>0
  row[trait+'_beta_alt']=-beta;row[trait+'_se']=se;row[trait+'_conditional_pip']=float(x[trait+'_conditional_abf_pip'])
 rows.append(row)
with (p/'alignment.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
(p/'receipt.json').write_text(json.dumps({'n':len(rows),'input_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (evidence,reporter)},'qualification':'Post-selected and conditional on reporter and source-allele conventions; no causal inference'},indent=2)+'\n');print(json.dumps(rows,indent=2))
