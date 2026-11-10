import csv,json,hashlib
from pathlib import Path
O=Path('outputs'); p1=O/'stage6_aou_eur_candidate_association.csv'; p2=O/'stage6_aou_eur_fatigue_candidate_association.csv'
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
ins={x['candidate']:x for x in read(p1) if x['phenotype']=='NS_333.2'}; fat={x['candidate']:x for x in read(p2)}
assert len(ins)==len(fat)==33 and set(ins)==set(fat)
out=[]
for v in sorted(ins):
 a,b=ins[v],fat[v]
 assert all(a[k]==b[k] for k in ('chrom','position','ref','alt','ancestry','sequencing_type','storage_source','status'))
 assert a['status']=='matched' and a['ancestry']=='eur' and b['phenotype']=='SS_807'
 x,y=float(a['pvalue']),float(b['pvalue']); assert 0<=x<=1 and 0<=y<=1
 out.append(dict(variant=v,chrom=a['chrom'],position38=a['position'],ref=a['ref'],alt=a['alt'],p_insomnia=x,p_fatigue=y,p_conjunction_max=max(x,y)))
indices=sorted(range(len(out)),key=lambda i:out[i]['p_conjunction_max']);q=1
for rank in range(len(out),0,-1):
 i=indices[rank-1]; q=min(q,out[i]['p_conjunction_max']*len(out)/rank);out[i]['bh_q_within_33']=q
p=O/'stage6_aou_insomnia_fatigue_eur_conjunction.csv'
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
summary={'n':len(out),'both_p_lt_05':sum(r['p_conjunction_max']<.05 for r in out),'bh_q_lt_05':sum(r['bh_q_within_33']<.05 for r in out),'smallest_max_p':min(r['p_conjunction_max'] for r in out),'best':min(out,key=lambda r:r['p_conjunction_max'])['variant']}
receipt={'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (p1,p2)},'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'summary':summary,'method':'Two-sided max(p_insomnia,p_fatigue), BH within 33 selected candidates','limits':'Postselected source and candidates, phenotype overlap, unknown source sample overlap, ancestry/reference and exact variant issues; not causal sharing or score validation.'}
(O/'stage6_aou_insomnia_fatigue_eur_conjunction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summary,indent=2))
