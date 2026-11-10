import csv,json,hashlib
from pathlib import Path
from statsmodels.stats.multitest import multipletests
O=Path('outputs');out=[];sources={}
for pop in ('eur','afr'):
 p=O/f'stage6_aou_{pop}_candidate_association.csv';sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();raw=list(csv.DictReader(p.open()));assert len(raw)==66
 groups={}
 for x in raw:
  assert x['status']=='matched' and x['ancestry']==pop and x['sequencing_type']=='genome' and x['storage_source']=='hail_gcs'
  d=groups.setdefault(x['candidate'],{});assert x['phenotype'] not in d;d[x['phenotype']]=x
 assert len(groups)==33
 these=[]
 for variant,d in groups.items():
  a=d['NS_333.2'];b=d['SS_807.21'];assert all(a[k]==b[k] for k in ('chrom','position','ref','alt'))
  p1,p2=float(a['pvalue']),float(b['pvalue']);assert 0<=p1<=1 and 0<=p2<=1
  these.append({'population':pop.upper(),'variant':variant,'chrom':a['chrom'],'position38':a['position'],'ref':a['ref'],'alt':a['alt'],'p_insomnia':p1,'p_weakness':p2,'p_conjunction_max':max(p1,p2)})
 for x,q in zip(these,multipletests([x['p_conjunction_max'] for x in these],method='fdr_bh')[1]):x['bh_q_within_33']=float(q)
 out.extend(these)
p=O/'stage6_aou_insomnia_weakness_conjunction.csv'
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
summary={}
for pop in ('EUR','AFR'):
 rows=[x for x in out if x['population']==pop];best=min(rows,key=lambda x:x['p_conjunction_max'])
 summary[pop]={'n':len(rows),'nominal_p_lt_05':sum(x['p_conjunction_max']<.05 for x in rows),'bh_q_lt_05':sum(x['bh_q_within_33']<.05 for x in rows),'smallest_p':best['p_conjunction_max'],'smallest_p_variant':best['variant']}
receipt={'source_sha256':sources,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'summary':summary,'method':'Intersection-union test max(p_insomnia,p_weakness); BH across 33 variants per group','limits':'Selected variants/outcomes and post hoc analysis; linkage disequilibrium and multiple comparisons beyond this family; uncertain source overlap; association with both traits does not prove shared causal variant, effect direction, score performance or frailty.'}
(O/'stage6_aou_insomnia_weakness_conjunction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
