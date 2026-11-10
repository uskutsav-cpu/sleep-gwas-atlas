#!/usr/bin/env python3
"""Join exact candidate associations from BBJ insomnia and Japanese knee GWAS.

This is a descriptive postselected join, not a regression or validation.
"""
import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
def read(name):
 with (p/name).open(newline='') as f:return list(csv.DictReader(f))
bfile='stage6_bbj_insomnia_oriented_candidates.csv';kfile='stage6_japanese_knee_candidate_matches.csv'
b={(x['score'],x['variant']):x for x in read(bfile)};k=read(kfile); out=[]
tfile='stage6_tpmi_shared22_insomnia_inventory.csv'
taiwan={x['variant']:x for x in read(tfile) if x['status']=='browser_lookup_rounded'}
for row in k:
 key=(row['score'],row['variant'])
 if row['score'] not in ('prototype15','shared22') or row['exact_rows']!='1':continue
 x=b[key]
 if x['status']!='exact position and pair':continue
 if (x['counted_allele'],x['other_allele'])!=(row['counted_allele'],row['other_allele']):raise ValueError('candidate allele mismatch '+repr(key))
 ib=float(x['beta_counted_assuming_saige']);ise=float(x['se']);kb=-float(row['beta_counted_if_source_coding_valid']);kse=float(row['source_standard_error'])
 t=taiwan.get(row['variant']) if row['score']=='shared22' else None
 if t and t['queue_identity_match']!='True':raise ValueError('Taiwan build/queue match missing '+repr(key))
 if t and ((t['counted_allele'],t['other_allele'])!=(row['counted_allele'],row['other_allele']) or {t['ref'],t['alt']}!={row['counted_allele'],row['other_allele']}):raise ValueError('Taiwan allele mismatch '+repr(key))
 tb=(float(t['insomnia_effect'])*(1 if t['counted_allele']==t['alt'] else -1)) if t else None
 out.append({'taiwan_insomnia_conditional_counted_beta':tb,'taiwan_insomnia_se_rounded':t['insomnia_se'] if t else None,'taiwan_insomnia_source_p_rounded':t['insomnia_p'] if t else None,'score':row['score'],'variant':row['variant'],'counted_allele':row['counted_allele'],'other_allele':row['other_allele'],'bbj_insomnia_beta':ib,'bbj_insomnia_se':ise,'bbj_insomnia_source_spa_p':x['source_p_spa'],'bbj_counted_allele_frequency':x['counted_allele_frequency'],'bbj_minor_dosage_count':x['minor_dosage_count'],'bbj_imputation_info':x['imputation_info'],'knee_lower_strength_beta':kb,'knee_standard_error':kse,'knee_source_p':row['source_p_value'],'knee_lower_strength_wald_low':kb-1.96*kse,'knee_lower_strength_wald_high':kb+1.96*kse,'direction_agreement':ib*kb>0})
assert len({(x['score'],x['variant']) for x in out})==len(out)
with (p/'stage6_bbj_knee_matched_candidates.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
counts={s:sum(x['score']==s for x in out) for s in ['prototype15','shared22']}
receipt={'inputs_sha256':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in [bfile,kfile,tfile]},'matched_counts':counts,'matched_shared22_taiwan':sum(x['score']=='shared22' and x['taiwan_insomnia_conditional_counted_beta'] is not None for x in out),'coding':'BBJ beta assumes Allele2 effect; knee higher-strength counted-allele beta negated; Taiwan beta assumes browser ALT is effect and is rounded','interpretation':'descriptive postselected comparison of distinct outcomes, not a joint phenotype, score regression, formal heterogeneity test, or predictive validation'}
(p/'stage6_bbj_knee_matched_candidates.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(counts))
