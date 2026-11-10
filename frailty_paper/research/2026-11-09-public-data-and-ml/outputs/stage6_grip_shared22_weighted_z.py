#!/usr/bin/env python3
"""Full 22-count marginal-Z screen, not patient-level score validation."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np

o=Path(__file__).resolve().parent
def read(name):
    with (o/name).open(newline='') as f:return list(csv.DictReader(f))
score_file=o/'unvalidated_22_shared_direction_score.csv'
score=read(score_file.name);grip={r['variant']:r for r in read('stage6_grip_gwas_candidate_aligned.csv') if r['score']=='shared22'}
freq={r['variant']:r for r in read('stage6_shared22_cross_population_frequencies.csv')}
covdata=json.loads((o/'shared22_count_variance_diagnostic.json').read_text())
source=json.loads((o/'stage6_grip_gwas_extract_receipt.json').read_text())
sha=hashlib.sha256(score_file.read_bytes()).hexdigest()
assert source['source_checksum_matches'] and covdata['score_sha256']==sha and len(score)==len(grip)==22
ids=[r['variant'] for r in score]; assert len(set(ids))==22 and len(covdata['pairs'])==231
C=np.diag([float(freq[v]['EUR_dosage_variance']) for v in ids]);indices={v:i for i,v in enumerate(ids)}
for p in covdata['pairs']:
    i,j=indices[p['variant_1']],indices[p['variant_2']];assert i!=j
    if p['same_chromosome']: C[i,j]=C[j,i]=float(p['dosage_covariance'])
sd_ref=np.sqrt(np.diag(C));assert np.all(sd_ref>0)
R=C/np.outer(sd_ref,sd_ref);assert np.linalg.eigvalsh(R).min()>-1e-8
z=[];a=[]
for s in score:
    g=grip[s['variant']];assert int(g['exact_rows'])==1 and s['counted_allele']==g['counted_allele'] and s['other_allele']==g['other_allele']
    p=float(g['counted_allele_frequency']);se=float(g['se']);assert 0<p<1 and se>0
    z.append(float(g['beta_counted'])/se);a.append(math.sqrt(2*p*(1-p)))
z=np.array(z);a=np.array(a)
rows=[]
for idx in [-1,*range(len(ids))]:
    keep=np.arange(len(ids))!=idx;w=a[keep];r=R[np.ix_(keep,keep)]
    stat=float(w@z[keep]/np.sqrt(w@r@w))
    rows.append({'omitted_variant':'none' if idx<0 else ids[idx],'n':int(sum(keep)),'weighted_z':stat,'nominal_two_sided_p':math.erfc(abs(stat)/math.sqrt(2))})
with (o/'stage6_grip_shared22_weighted_z.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
receipt={'score_sha256':sha,'gwas_md5':source['observed_md5'],'n_target_variants':22,'reference':'1000 Genomes phase 3 EUR503 dosage covariance; cross-chromosome values set to zero','statistic':'a.T z / sqrt(a.T R a); a=sqrt(2*target_EAF*(1-target_EAF))','results':rows[0],'leave_one_out_min':min(x['weighted_z'] for x in rows[1:]),'leave_one_out_max':max(x['weighted_z'] for x in rows[1:]),'limitations':'Postselected; target contains UK Biobank discovery participants; mixed-cohort meta-analysis with varying variant-level participation; no independent prediction, discrimination or calibration.'}
(o/'stage6_grip_shared22_weighted_z.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
