import csv,json,hashlib
from pathlib import Path
from statsmodels.stats.multitest import multipletests
p=Path('outputs/stage6_aou_full_candidate_association.csv')
r=list(csv.DictReader(p.open()));assert len(r)==66 and all(x['status']=='matched' for x in r)
for phenotype in ('NS_333.2','SS_807.21'):
 rows=[x for x in r if x['phenotype']==phenotype];assert len(rows)==33
 for x,q in zip(rows,multipletests([float(v['pvalue']) for v in rows],method='fdr_bh')[1]):x['bh_q_33_within_phenotype']=float(q)
for x,q in zip(r,multipletests([float(x['pvalue']) for x in r],method='fdr_bh')[1]):x['bh_q_66_joint']=float(q)
out=Path('outputs/stage6_aou_full_candidate_association_with_fdr.csv')
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)
summary={}
for ph in ('NS_333.2','SS_807.21'):
 rows=[x for x in r if x['phenotype']==ph]
 summary[ph]={'n':len(rows),'nominal_p_lt_05':sum(float(x['pvalue'])<.05 for x in rows),'bh_q_33_lt_05':sum(x['bh_q_33_within_phenotype']<.05 for x in rows),'bh_q_66_lt_05':sum(x['bh_q_66_joint']<.05 for x in rows),'bh_33_hits':[{'variant':x['candidate'],'p':float(x['pvalue']),'q':x['bh_q_33_within_phenotype']} for x in rows if x['bh_q_33_within_phenotype']<.05]}
res={'summary':summary,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'limits':'Descriptive, post-selected variants and outcomes; BH family excludes broader project exploration, and association p-values do not validate scores.'}
Path('outputs/stage6_aou_full_candidate_association_fdr_receipt.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
