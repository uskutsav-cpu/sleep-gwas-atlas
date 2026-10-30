#!/usr/bin/env python3
"""Independent rowwise reconstruction and original-file checks, separate from baseline scoring."""
from pathlib import Path
import math
import numpy as np,pandas as pd,pyreadstat
P=Path(__file__).resolve().parent;R=P/'source_files'
x=pd.read_csv(P/'finalized_participant_data.csv');raw=pd.read_csv(P/'raw_variable_dictionary.csv');D=pd.read_csv(P/'fi_item_dictionary.csv');assert len(D)==49
checks=[]
# Verify all original raw values on age 60+ by rereading each CDC source, including special missing codes.
for (year,f),group in raw.groupby(['cycle','source_file']):
 a,_=pyreadstat.read_xport(str(R/f),encoding='latin1');names=group.variable.tolist();names=[n for n in names if n in x and n!='SEQN']
 if f.startswith('RXQ_RX'):a=a.drop_duplicates('SEQN')
 assert a.SEQN.is_unique
 cycle=f'{year}-{year+1}';q=x[x.cycle==cycle][['SEQN']+names].merge(a[['SEQN']+names],on='SEQN',how='left',validate='one_to_one',suffixes=('_clean','_orig'))
 for v in names:
  a1=q[v+'_clean'];a2=q[v+'_orig'];bad=~(a1.eq(a2)|(a1.isna()&a2.isna()));checks.append({'kind':'raw_source','item':v,'cycle':cycle,'mismatches':int(bad.sum()),'compared':len(q)})

def yes(v):return 1. if v==1 else 0. if v==2 else np.nan
functions={'money':'A','stoop':'D','lift':'E','chores':'F','meals':'G','chair':'I','bed':'J','eat':'K','dress':'L','stand':'M','grasp':'P','social':'R','push':'T','walk':'B','steps':'C'}
mood={'interest':'DPQ010','depressed':'DPQ020','sleep':'DPQ030','fatigue':'DPQ040','appetite':'DPQ050','self_worth':'DPQ060','concentration':'DPQ070'}
disease={'arthritis':'MCQ160A','thyroid':'MCQ160M','bronchitis':'MCQ160K','cancer':'MCQ220','heart_failure':'MCQ160B','coronary':'MCQ160C','angina':'MCQ160D','heart_attack':'MCQ160E','stroke':'MCQ160F','hypertension':'BPQ020','kidney':'KIQ022'}
def rec(row):
 out={'memory':yes(row.PFQ057)}
 for nm,code in functions.items():
  q=row['PFQ061'+code];out[nm]=(0. if q==1 else 1. if q in (2,3,4) else np.nan)
  if code in ['B','C'] and row.PFQ054==1:out[nm]=1.
 for nm,v in mood.items():
  q=row[v];out[nm]=q/3 if q in (0,1,2,3) else np.nan
 for nm,v in disease.items():out[nm]=yes(row[v])
 out['diabetes']={1:1.,2:0.,3:.5}.get(row.DIQ010,np.nan);out['urine_leak']=0. if row.KIQ005==1 else 1. if row.KIQ005 in (2,3,4,5) else np.nan
 out['general_health']=0. if row.HUQ010 in (1,2,3) else 1. if row.HUQ010 in (4,5) else np.nan
 out['health_worse']=1. if row.HUQ020==2 else 0. if row.HUQ020 in (1,3) else np.nan
 out['hospital']=yes(row.HUQ071)
 if row.cycle=='2011-2012': out['care_visits']={0:0.,1:.5,2:.5,4:1.,5:1.}.get(row.HUQ050,np.nan)
 else:out['care_visits']={0:0.,1:.5,2:.5,3:.5,4:1.,5:1.,6:1.,7:1.,8:1.}.get(row.HUQ051,np.nan)
 out['medicines']=0. if row.RXDUSE==2 else (.5 if 1<=row.RXDCOUNT<=4 else 1. if 5<=row.RXDCOUNT<=23 else np.nan) if row.RXDUSE==1 else np.nan
 bmi=row.BMXBMI;out['bmi']=np.nan if pd.isna(bmi) else 1. if bmi<18.5 or bmi>=30 else .5 if bmi>=25 else 0.
 grip=np.nan
 if row.MGDEXSTS==1 and row.MGD130 in (1,2) and row.MGATHAND in (1,2):
  hand=1 if row.MGD130==row.MGATHAND else 2;trials=[row[f'MGXH{hand}T{i}'] for i in (1,2,3)];valid=[v for v in trials if pd.notna(v)];grip=sum(valid)/len(valid) if valid else np.nan
 if pd.isna(grip) or pd.isna(bmi) or row.RIAGENDR not in (1,2):out['weak_grip']=np.nan
 else:
  if row.RIAGENDR==1:cut=29 if bmi<=24 else 30 if bmi<=28 else 32
  else:cut=17 if bmi<=23 else 17.3 if bmi<=26 else 18 if bmi<=29 else 21
  out['weak_grip']=float(grip<=cut)
 for nm,v in [('hba1c','LBXGH'),('rbc','LBXRBCSI'),('hemoglobin','LBXHGB'),('rdw','LBXRDW'),('lymphocytes','LBXLYPCT'),('neutrophils','LBXNEPCT')]:
  q=row[v]
  if pd.isna(q):out[nm]=np.nan;continue
  ranges={'rbc':(4.7,6.1) if row.RIAGENDR==1 else (4.2,5.4),'hemoglobin':(13.5,18) if row.RIAGENDR==1 else (12,16),'rdw':(11.6,14.6),'lymphocytes':(20,40),'neutrophils':(40,80)}
  out[nm]=float(q>5.7) if nm=='hba1c' else float(not(ranges[nm][0]<=q<=ranges[nm][1]))
 return out
rep=pd.DataFrame([rec(row) for _,row in x.iterrows()]).add_prefix('fi_')
for v in D.item:
 a=x[v].to_numpy();b=rep[v].to_numpy();bad=~np.isclose(a,b,atol=1e-12,rtol=0,equal_nan=True);checks.append({'kind':'independent_recode','item':v,'cycle':'pooled','mismatches':int(bad.sum()),'compared':len(x)})
# Independent exclusion: exact sets and arithmetic for each version; no raw directly sleep/fatigue question appears in neither-index.
items=D.item.tolist();versions={'fi_all':items,'fi_no_sleep':[i for i in items if i!='fi_sleep'],'fi_no_fatigue':[i for i in items if i!='fi_fatigue'],'fi_neither':[i for i in items if i not in ['fi_sleep','fi_fatigue']]}
for name,cols in versions.items():
 assert ('fi_sleep' in cols)==(name in ('fi_all','fi_no_fatigue'));assert ('fi_fatigue' in cols)==(name in ('fi_all','fi_no_sleep'))
 n=rep[cols].notna().sum(axis=1);score=(rep[cols].sum(axis=1,min_count=1)/n).where(n>=math.ceil(.8*len(cols)));bad=~np.isclose(score,x[name],atol=1e-12,rtol=0,equal_nan=True);checks.append({'kind':'independent_score','item':name,'cycle':'pooled','mismatches':int(bad.sum()),'compared':len(x)})
pd.DataFrame(checks).to_csv(P/'supplement_item_audit_checks.csv',index=False);err=sum(r['mismatches'] for r in checks);print('audit checks',len(checks),'mismatches',err);assert err==0
