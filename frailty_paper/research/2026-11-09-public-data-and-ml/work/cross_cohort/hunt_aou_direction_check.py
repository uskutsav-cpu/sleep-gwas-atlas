#!/usr/bin/env python3
import csv,hashlib,json,re
from pathlib import Path
root=Path(__file__).resolve().parents[2]
p1=root/'outputs/stage6_lane_aou_six_variant_comparison.csv';p2=root/'outputs/stage6_aou_eur_candidate_association.csv'
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
aou={(r['candidate'],r['phenotype']):r for r in read(p2)};out=[]
for r in read(p1):
 x=aou[r['variant'],'NS_333.2'];e,o=r['lane_effect_other'].split('/');ref,alt=x['ref'],x['alt']
 if {e,o}!={ref,alt}:raise ValueError(r['variant'])
 beta=float(x['beta']);or_hunt=float(re.match(r'^([0-9.]+)',r['hunt_or_ci']).group(1))
 aligned=beta if e==alt else -beta
 out.append({'variant':r['variant'],'hunt_effect_other':r['lane_effect_other'],'hunt_or_ci':r['hunt_or_ci'],'hunt_p':r['hunt_p'],'aou_ref_alt':ref+'/'+alt,'aou_beta_raw':beta,'aou_beta_on_hunt_effect_assuming_alt':aligned,'aou_p':x['pvalue'],'direction_agrees_assuming_alt':str((or_hunt-1)*aligned>0).lower(),'caveat':'Conditional allele-coding assumption; source OR is rounded; variant identity, sample overlap and phenotype matching remain unverified.'})
outp=root/'outputs/stage6_hunt_aou_allele_direction_sensitivity.csv'
with outp.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print('rows',len(out),'agreement',sum(r['direction_agrees_assuming_alt']=='true' for r in out));print('output_sha256',hashlib.sha256(outp.read_bytes()).hexdigest())
(root/'outputs/stage6_hunt_aou_allele_direction_receipt.json').write_text(json.dumps({'inputs':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [p1,p2]},'output_sha256':hashlib.sha256(outp.read_bytes()).hexdigest(),'rows':len(out)},indent=2)+'\n')
