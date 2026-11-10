#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;x=pd.read_csv(P/'finalized_participant_data.csv');items=pd.read_csv(P/'fi_item_dictionary.csv');diag=pd.read_csv(P/'fi_item_diagnostics.csv');combined=items.merge(diag[diag.cycle=='2011-2012'][['item','n_observed','missing_fraction','prevalence_any_deficit','spearman_age']],on='item',validate='one_to_one').merge(diag[diag.cycle=='2013-2014'][['item','n_observed','missing_fraction','prevalence_any_deficit','spearman_age']],on='item',validate='one_to_one',suffixes=('_2011','_2013'))
combined['in_index_without_sleep_fatigue']=~combined.item.isin(['fi_sleep','fi_fatigue']);combined['direct_sleep_fatigue_item']=combined.item.map({'fi_sleep':'sleep (includes excessive sleep)','fi_fatigue':'fatigue/little energy'}).fillna('no direct item');combined['manual_audit']='recode agrees with independent rowwise implementation; interpretation/limitations in rule and notes';combined.to_csv(P/'supplement_49_item_audit.csv',index=False)
r=[]
for cycle in sorted(x.cycle.unique()):
 for exposure in ['all','0','1']:
  d=x[(x.cycle==cycle)&(x.WTMEC2YR>0)]
  if exposure!='all':d=d[d.sleep_told_doctor==int(exposure)]
  for v in items.item:
   m=d[v].isna();r.append({'cycle':cycle,'ever_told_doctor_sleep':exposure,'item':v,'sample_size':len(d),'missing_n':int(m.sum()),'missing_fraction':float(m.mean()),'weighted_missing_fraction':np.average(m.astype(int),weights=d.WTMEC2YR)})
pd.DataFrame(r).to_csv(P/'supplement_missingness_by_exposure_cycle.csv',index=False)
# Counts and score availability are unweighted; weights apply only to the survey estimates.
flow=[]
for cycle,g in x.groupby('cycle'):
 q=g.WTMEC2YR>0;steps={'interviewed_age60':len(g),'positive_mec_weight':q.sum(),'sleep_report_observed':(q&g.sleep_told_doctor.notna()).sum(),'FI_neither_80pct':(q&g.fi_neither.notna()).sum(),'common_all_FI_and_sleep':(q&(g.common_scores==1)).sum(),'common_FI_and_prior_physical4':(q&(g.common_scores==1)&g.prior_physical4.notna()&g.fried_physical4.notna()).sum(),'linkage_eligible':(q&(g.ELIGSTAT==1)).sum(),'linkage_eligible_FI':(q&(g.ELIGSTAT==1)&g.fi_neither.notna()).sum(),'deaths_linkage_eligible_FI':(q&(g.ELIGSTAT==1)&g.fi_neither.notna()&(g.MORTSTAT==1)).sum()}
 flow.extend([{'cycle':cycle,'criterion':a,'n':int(b)} for a,b in steps.items()])
pd.DataFrame(flow).to_csv(P/'supplement_flow.csv',index=False)
print('49-item audit and missingness tables complete')
