#!/usr/bin/env python3
"""Read the preserved reconstruction and CDC public source files; write new research datasets."""
from pathlib import Path
import hashlib, numpy as np, pandas as pd, pyreadstat
P=Path(__file__).resolve().parent;R=P/'source_files'
# baseline_analysis.py writes this local copy, never the earlier project.
x=pd.read_csv(P/'nhanes_2011_2014_clean.csv');dict=pd.read_csv(P/'fi_item_dictionary.csv');items=dict.item.tolist();assert len(items)==49 and len(set(items))==49
systems={'sleepfatigue':{'fi_sleep','fi_fatigue'},'mood':set(dict.loc[dict.domain=='mood','item']),'lab':set(dict.loc[dict.domain=='laboratory','item']),'care':set(dict.loc[dict.domain=='healthcare','item'])}
versions={'fi_no_mood_no_care':[v for v in items if v not in systems['mood']|systems['care']], 'fi_neither_no_labs':[v for v in items if v not in systems['sleepfatigue']|systems['lab']], 'fi_neither_no_function_walk':[v for v in items if v not in systems['sleepfatigue']|{'fi_walk','fi_steps'}], 'fi_neither_no_bmi':[v for v in items if v not in systems['sleepfatigue']|{'fi_bmi'}]}
for name,cols in versions.items():
 n=x[cols].notna().sum(axis=1);x[name+'_n']=n;x[name]=(x[cols].sum(axis=1,min_count=1)/n).where(n>=np.ceil(.8*len(cols)))
# Alternative exposure uses DPQ030, which is categorically excluded from the outcome.
x['sleep_recent_halfdays']=np.where(x.DPQ030.isin([0,1,2,3]),np.where(x.DPQ030>=2,1,0),np.nan)
x['sleep_recent_any']=np.where(x.DPQ030.isin([0,1,2,3]),np.where(x.DPQ030>=1,1,0),np.nan)
x['FI_eligible']=x.fi_neither.notna().astype(int);x['common_scores']=x[['fi_all','fi_no_sleep','fi_no_fatigue','fi_neither','sleep_told_doctor']].notna().all(axis=1).astype(int)
x['fi_difference']=x.fi_all-x.fi_neither;x['fi_high']=np.where(x.fi_neither.notna(),(x.fi_neither>.21).astype(int),np.nan)
# Public mortality record, eligibility and vital status only: no CDC-perturbed times are analyzed.
morts=[];full=[]
for y,s in [(2011,'G'),(2013,'H')]:
 cycle=f'{y}-{y+1}';f=R/f'NHANES_{y}_{y+1}_MORT_2019_PUBLIC.dat';
 if not f.exists():raise FileNotFoundError(f)
 m=pd.read_fwf(f,colspecs=[(0,6),(14,15),(15,16)],names=['SEQN','ELIGSTAT','MORTSTAT'],na_values=['','.']);m['cycle']=cycle;assert not m.SEQN.duplicated().any();morts.append(m)
 d,_=pyreadstat.read_xport(str(R/f'DEMO_{s}.XPT'),encoding='latin1');d=d[['SEQN','RIDAGEYR','WTMEC2YR','SDMVSTRA','SDMVPSU']];d['cycle']=cycle;full.append(d)
m=pd.concat(morts,ignore_index=True);x=x.merge(m,on=['cycle','SEQN'],how='left',validate='one_to_one'); assert x.ELIGSTAT.notna().all()
x.to_csv(P/'finalized_participant_data.csv',index=False)
D=pd.concat(full,ignore_index=True);D=D[(D.WTMEC2YR>0)&D.SDMVSTRA.notna()&D.SDMVPSU.notna()].copy(); cols=[v for v in x if v not in ['RIDAGEYR','WTMEC2YR','SDMVSTRA','SDMVPSU']];D=D.merge(x[cols],on=['cycle','SEQN'],how='left',validate='one_to_one');D['wt4']=D.WTMEC2YR/2;D['stratum']=D.cycle+':'+D.SDMVSTRA.astype(int).astype(str);D.to_csv(P/'survey_full_design_frame.csv',index=False)
pd.DataFrame([{'score':n,'n_items':len(c),'items':';'.join(c),'minimum_observed':int(np.ceil(.8*len(c)))} for n,c in versions.items()]).to_csv(P/'alternative_score_definitions.csv',index=False)
pd.DataFrame([{'file':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'url':'https://ftp.cdc.gov/pub/health_statistics/NCHS/datalinkage/linked_mortality/'+f.name} for f in sorted(R.glob('*MORT_2019_PUBLIC.dat'))]).to_csv(P/'mortality_source_checksums.csv',index=False)
print('older',len(x),'MEC',int((x.WTMEC2YR>0).sum()),'FI common',int(((x.WTMEC2YR>0)&(x.common_scores==1)).sum()),'linked',int(((x.WTMEC2YR>0)&(x.ELIGSTAT==1)).sum()))
