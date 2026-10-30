#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent
x=pd.read_csv(P/'finalized_participant_data.csv');r=pd.read_csv(P/'supplement_survey_models.csv');prior=pd.read_csv(P/'survey_results.csv')
for new,old in [('fi_all','fi_all'),('fi_no_sleep','fi_no_sleep'),('fi_no_fatigue','fi_no_fatigue'),('fi_neither','fi_neither'),('fi_difference','paired_fi_all_minus_fi_neither'),('fi_diff_nosleep','paired_fi_all_minus_fi_no_sleep'),('fi_diff_nofatigue','paired_fi_all_minus_fi_no_fatigue')]:
 a=r[(r.analysis==new)&(r.group=='pooled')].iloc[0];b=prior[(prior.outcome==old)&(prior.adjustment=='age_sex')&(prior.cycle=='pooled')].iloc[0];assert a.n==b.n and abs(a.estimate-b.estimate)<1e-10 and abs(a.se_link-b.se)<1e-10,(new,a.estimate,b.estimate,a.se_link,b.se)
assert len(r[r.family=='primary'])==7 and len(r[r.family=='exploratory'])==33
assert pd.read_csv(P/'supplement_item_audit_checks.csv').mismatches.sum()==0
for n in (15,22):
 s=pd.read_csv(P/f'score_{n}_variants_standardized.csv');assert len(s)==n and s.rsid.nunique()==n
k=pd.read_csv(P/'supplement_genetic_reference_audit.csv');assert len(k)==37 and k.rsid.nunique()==33 and (k.request_status==200).all() and (k.exact_locus_mappings==1).all() and k.score_pair_present.all()
w=pd.read_csv(P/'supplement_gwas_weight_audit.csv');assert len(w)==111 and (w.status=='checked').all() and w.weight_agrees.all()
for cyc,part in x.groupby('cycle'):
 name='NHANES_'+cyc.replace('-','_')+'_MORT_2019_PUBLIC.dat';a=pd.read_fwf(P/'source_files'/name,colspecs=[(0,6),(14,15),(15,16)],names=['SEQN','e','m'],na_values=['','.']);j=part.merge(a,on='SEQN',how='left',validate='one_to_one');assert (j.ELIGSTAT.eq(j.e)|(j.ELIGSTAT.isna()&j.e.isna())).all();assert (j.MORTSTAT.eq(j.m)|(j.MORTSTAT.isna()&j.m.isna())).all()
for _,s in pd.read_csv(P/'source_checksums.csv').iterrows():assert hashlib.sha256((P/'source_files'/s.file).read_bytes()).hexdigest()==s.sha256
for _,s in pd.read_csv(P/'mortality_source_checksums.csv').iterrows():assert hashlib.sha256((P/'source_files'/s.file).read_bytes()).hexdigest()==s.sha256
for n in ['figure1_fi_associations','figure2_distribution_mortality']:
 for ext in ['png','svg','pdf']:assert (P/f'{n}.{ext}').stat().st_size>1000
print('PASS: independent 49-item and original-source audit; all 7 R/Python survey models; unperturbed mortality record linkage; 37 reference/111 summary-weight checks; source hashes; figure formats')
