import csv,json,hashlib
from pathlib import Path
from statsmodels.stats.multitest import multipletests
pairs=[('outputs/stage6_aou_eur_walking_candidate_association.csv','pvalue','outputs/stage6_aou_eur_walking_candidate_fdr.csv',33),('outputs/stage6_aou_eur_walking_weighted_z_sensitivity.csv','nominal_p_two_sided','outputs/stage6_aou_eur_walking_score_fdr.csv',7)]
receipt={}
for src,col,dest,n in pairs:
 rows=list(csv.DictReader(open(src)));assert len(rows)==n
 ps=[float(r[col]) for r in rows];assert all(0<=p<=1 for p in ps)
 for r,q in zip(rows,multipletests(ps,method='fdr_bh')[1]):r['bh_q_within_family']=float(q)
 with open(dest,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 receipt[dest]={'rows':n,'nominal_p_lt_05':sum(p<.05 for p in ps),'bh_q_lt_05':sum(float(r['bh_q_within_family'])<.05 for r in rows),'source_sha256':hashlib.sha256(Path(src).read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(Path(dest).read_bytes()).hexdigest()}
Path('outputs/stage6_aou_eur_walking_fdr_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
