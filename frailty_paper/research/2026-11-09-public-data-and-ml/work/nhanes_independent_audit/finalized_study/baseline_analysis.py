#!/usr/bin/env python3
"""Reconstruct an adapted 49-deficit NHANES FI and Fried-like components, 2011–2014."""
from pathlib import Path
from urllib.request import urlopen
import hashlib, ssl, certifi
import numpy as np, pandas as pd, pyreadstat
from scipy import stats
from bs4 import BeautifulSoup
P=Path(__file__).resolve().parent; R=P/'source_files';R.mkdir(exist_ok=True)
F={
'DEMO':'SEQN RIDAGEYR RIAGENDR RIDRETH3 DMDEDUC2 INDFMPIR WTMEC2YR SDMVSTRA SDMVPSU RIDSTATR',
'SLQ':'SEQN SLQ050', 'DPQ':'SEQN DPQ010 DPQ020 DPQ030 DPQ040 DPQ050 DPQ060 DPQ070',
'PFQ':'SEQN PFQ054 PFQ057 PFQ061A PFQ061B PFQ061C PFQ061D PFQ061E PFQ061F PFQ061G PFQ061I PFQ061J PFQ061K PFQ061L PFQ061M PFQ061P PFQ061R PFQ061T',
'MCQ':'SEQN MCQ160A MCQ160B MCQ160C MCQ160D MCQ160E MCQ160F MCQ160K MCQ160M MCQ220',
'BPQ':'SEQN BPQ020','DIQ':'SEQN DIQ010','KIQ_U':'SEQN KIQ022 KIQ005','HUQ':'SEQN HUQ010 HUQ020 HUQ071','RXQ_RX':'SEQN RXDUSE RXDCOUNT',
'CBC':'SEQN LBXLYPCT LBXNEPCT LBXRBCSI LBXHGB LBXRDW','GHB':'SEQN LBXGH',
'MGX':'SEQN MGD130 MGATHAND MGDEXSTS MGDCGSZ MGXH1T1 MGXH1T2 MGXH1T3 MGXH2T1 MGXH2T2 MGXH2T3 MGXH1T1E MGXH1T2E MGXH1T3E MGXH2T1E MGXH2T2E MGXH2T3E',
'BMX':'SEQN BMXBMI BMXWT BMXHT', 'WHQ':'SEQN WHD020 WHD050 WHQ060 WHQ070',
'PAQ':'SEQN PAQ605 PAQ610 PAD615 PAQ620 PAQ625 PAD630 PAQ635 PAQ640 PAD645 PAQ650 PAQ655 PAD660 PAQ665 PAQ670 PAD675', 'SMQ':'SEQN SMQ020'}
rawdict=[]; checks=[]; frames=[]
for yr,suffix in [(2011,'G'),(2013,'H')]:
 base=f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr}/DataFiles/'; dat=None
 for f,columns in F.items():
  for ext in ['XPT','htm']:
   p=R/f'{f}_{suffix}.{ext}';url=base+p.name
   if not p.exists():
    with urlopen(url,context=ssl.create_default_context(cafile=certifi.where())) as resp:p.write_bytes(resp.read())
   checks.append({'cycle':str(yr)+'-'+str(yr+1),'file':p.name,'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
  frame,meta=pyreadstat.read_xport(str(R/f'{f}_{suffix}.XPT'),encoding='latin1'); columns=columns.split(); columns += ([('HUQ050' if yr==2011 else 'HUQ051')] if f=='HUQ' else []);assert set(columns)<=set(frame.columns),(yr,f,set(columns)-set(frame.columns))
  for c in columns:
   if c=='SEQN' and f!='DEMO':continue
   rawdict.append({'cycle':yr,'variable':c,'source_file':f'{f}_{suffix}.XPT','description':meta.column_names_to_labels.get(c,''),'documentation_url':base+f'{f}_{suffix}.htm','raw_coding':'Unmodified CDC values; see codebook for missing and special codes'})
  frame=frame[columns]
  if f=='RXQ_RX':
   assert (frame.groupby('SEQN')[['RXDUSE','RXDCOUNT']].nunique(dropna=True)<=1).all().all(),'Inconsistent participant prescription count'
   frame=frame.drop_duplicates('SEQN')
  assert frame.SEQN.is_unique,(yr,f)
  dat=frame if dat is None else dat.merge(frame,on='SEQN',how='left',validate='one_to_one')
 dat['cycle']=f'{yr}-{yr+1}';frames.append(dat)
d=pd.concat(frames,ignore_index=True);assert not d.duplicated(['cycle','SEQN']).any()
pd.DataFrame(checks).to_csv(P/'source_checksums.csv',index=False);pd.DataFrame(rawdict).to_csv(P/'raw_variable_dictionary.csv',index=False)
def numeric_valid(x,codes):return x.where(x.isin(codes))
def binary_yes(x):return numeric_valid(x,[1,2]).map({1:1.,2:0.})
def numcond(x,cond):return np.where(x.notna(),cond.astype(float),np.nan)
itemdict=[];items=[]
def item(name,value,domain,source,rule,notes=''):
 col='fi_'+name;d[col]=value;items.append(col);itemdict.append(dict(item=col,domain=domain,source=source,rule=rule,notes=notes))
item('memory',binary_yes(d.PFQ057),'cognition','PFQ057','1=yes deficit; 2=no; 7/9/blank missing')
for letter,label in [('A','money'),('D','stoop'),('E','lift'),('F','chores'),('G','meals'),('I','chair'),('J','bed'),('K','eat'),('L','dress'),('M','stand'),('P','grasp'),('R','social'),('T','push'),('B','walk'),('C','steps')]:
 v='PFQ061'+letter;x=d[v];a=np.where(x.isin([1,2,3,4]),(x>=2).astype(float),np.nan)
 if letter in 'BC':a=np.where(d.PFQ054==1,1,a)
 item(label,a,'function',v+(' PFQ054' if letter in 'BC' else ''),'1=no difficulty ->0; 2-4 ->1; 5/7/9/blank missing'+('; PFQ054=1 ->1 as proxy for skipped item' if letter in 'BC' else ''),'Quarter-mile and steps skipped if equipment needed; infer limitation, not actual task ability.' if letter in 'BC' else '')
for v,label in [('DPQ010','interest'),('DPQ020','depressed'),('DPQ030','sleep'),('DPQ040','fatigue'),('DPQ050','appetite'),('DPQ060','self_worth'),('DPQ070','concentration')]:item(label,numeric_valid(d[v],[0,1,2,3])/3,'mood',v,'0,1,2,3 -> 0, 1/3, 2/3, 1; 7/9/blank missing')
for v,label in [('MCQ160A','arthritis'),('MCQ160M','thyroid'),('MCQ160K','bronchitis'),('MCQ220','cancer'),('MCQ160B','heart_failure'),('MCQ160C','coronary'),('MCQ160D','angina'),('MCQ160E','heart_attack'),('MCQ160F','stroke'),('BPQ020','hypertension'),('KIQ022','kidney')]:item(label,binary_yes(d[v]),'comorbidity',v,'1=yes ->1; 2=no ->0; other codes missing')
item('diabetes',numeric_valid(d.DIQ010,[1,2,3]).map({1:1.,2:0.,3:.5}),'comorbidity','DIQ010','1=yes ->1; 2=no ->0; 3 borderline ->0.5; other missing')
item('urine_leak',numeric_valid(d.KIQ005,[1,2,3,4,5]).map({1:0.,2:1.,3:1.,4:1.,5:1.}),'comorbidity','KIQ005','1 never ->0; 2-5 any leakage ->1; other missing','Paper did not specify item or frequency threshold.')
item('general_health',np.where(d.HUQ010.isin([1,2,3,4,5]),(d.HUQ010>=4).astype(float),np.nan),'healthcare','HUQ010','1-3 ->0; 4-5 ->1; 7/9/blank missing')
item('health_worse',numeric_valid(d.HUQ020,[1,2,3]).map({1:0.,2:1.,3:0.}),'healthcare','HUQ020','2 worse ->1; 1 better or 3 same ->0; other missing')
item('hospital',binary_yes(d.HUQ071),'healthcare','HUQ071','1=yes ->1; 2=no ->0; other missing')
item('care_visits',np.where(d.cycle=='2011-2012',numeric_valid(d.HUQ050,[0,1,2,3,4,5]).map({0:0.,1:.5,2:.5,3:np.nan,4:1.,5:1.}),numeric_valid(d.HUQ051,list(range(9))).map({0:0.,1:.5,2:.5,3:.5,4:1.,5:1.,6:1.,7:1.,8:1.})),'healthcare','HUQ050 (2011) HUQ051 (2013)','0 visits ->0; 1-5 ->0.5; 6+ ->1; 2011 4-9 category unresolved and missing','2011 includes emergency/home visits; 2013 excludes them; 2011 4–9 spans threshold.')
med=np.where(d.RXDUSE==2,0,np.where((d.RXDUSE==1)&d.RXDCOUNT.between(1,23),np.where(d.RXDCOUNT>=5,1,.5),np.nan))
item('medicines',med,'healthcare','RXDUSE RXDCOUNT','RXDUSE=2 ->0; RXDUSE=1 and RXDCOUNT 1-4 ->0.5, 5+ ->1; otherwise missing')
bmi=d.BMXBMI;item('bmi',np.where(bmi.notna(),np.where((bmi<18.5)|(bmi>=30),1,np.where(bmi>=25,.5,0)),np.nan),'performance','BMXBMI','BMI <18.5 or >=30 ->1; 25 to <30 ->0.5; 18.5 to <25 ->0','Published table says <18.5, 30; interpreted >=30 as obesity.')
# NHANES starts with MGATHAND; H1 is first tested. Use mean of available dominant hand trials, both hands required to avoid ambiguous single-hand indexing.
h1=d[['MGXH1T1','MGXH1T2','MGXH1T3']].mean(axis=1);h2=d[['MGXH2T1','MGXH2T2','MGXH2T3']].mean(axis=1)
d['dominant_grip_mean']=np.where((d.MGDEXSTS==1)&d.MGD130.isin([1,2])&d.MGATHAND.isin([1,2]),np.where(d.MGD130==d.MGATHAND,h1,h2),np.nan)
male=d.RIAGENDR==1;female=d.RIAGENDR==2
cut=np.where(male,np.select([bmi<=24,bmi<=28],[29,30],default=32),np.select([bmi<=23,bmi<=26,bmi<=29],[17,17.3,18],default=21))
d['grip_threshold_kg']=np.where(bmi.notna()&(male|female),cut,np.nan)
item('weak_grip',np.where(d.dominant_grip_mean.notna()&d.grip_threshold_kg.notna(),(d.dominant_grip_mean<=d.grip_threshold_kg).astype(float),np.nan),'performance','MGD130 MGATHAND MGXH1T1-T3 MGXH2T1-T3 BMXBMI RIAGENDR','Mean available dominant-hand trials <= sex/BMI-specific Fried thresholds; both hands tested, known dominance required','NHANES dynamometer and procedure differ; mapping H1 to first tested hand is inferred from protocol; published FI does not specify hand/aggregation.')
for name,v,ok,rule in [('hba1c','LBXGH',(d.LBXGH<=5.7),'<=5.7 normal'),('rbc','LBXRBCSI',np.where(male,d.LBXRBCSI.between(4.7,6.1),d.LBXRBCSI.between(4.2,5.4)),'male 4.7-6.1; female 4.2-5.4'),('hemoglobin','LBXHGB',np.where(male,d.LBXHGB.between(13.5,18),d.LBXHGB.between(12,16)),'male 13.5-18; female 12-16'),('rdw','LBXRDW',d.LBXRDW.between(11.6,14.6),'11.6-14.6 normal'),('lymphocytes','LBXLYPCT',d.LBXLYPCT.between(20,40),'20-40 normal'),('neutrophils','LBXNEPCT',d.LBXNEPCT.between(40,80),'40-80 normal')]:item(name,np.where(d[v].notna(),1-np.asarray(ok,dtype=int),np.nan),'laboratory',v,rule+' ->0; outside ->1; missing ->missing','Published reference intervals are study conventions, not clinical diagnosis.')
assert len(items)==49
versions={'fi_all':items,'fi_no_sleep':[x for x in items if x!='fi_sleep'],'fi_no_fatigue':[x for x in items if x!='fi_fatigue'],'fi_neither':[x for x in items if x not in ['fi_sleep','fi_fatigue']],'fi_no_mood':[x for x in items if x not in ['fi_interest','fi_depressed','fi_sleep','fi_fatigue','fi_appetite','fi_self_worth','fi_concentration']],'fi_neither_no_care':[x for x in items if x not in ['fi_sleep','fi_fatigue','fi_care_visits']]}
# Secondary data-quality screen: <=5% missing in each older-adult cycle. This can yield fewer than 30 items.
ref=d.loc[(d.RIDAGEYR>=60)&(d.WTMEC2YR>0)]; strictitems=[v for v in items if all(ref.loc[ref.cycle==cy,v].isna().mean()<=.05 for cy in ['2011-2012','2013-2014'])]; versions['fi_missing5_subset']=strictitems
for n,cols in versions.items():
 d[n+'_n']=d[cols].notna().sum(axis=1);d[n]=d[cols].sum(axis=1,min_count=1)/d[n+'_n'];d.loc[d[n+'_n']<int(np.ceil(.8*len(cols))),n]=np.nan
 d[n+'_complete']=d[cols].mean(axis=1).where(d[n+'_n']==len(cols))
# Original Fried shrinking based on recalled loss >10 lb, and WHQ060 intent response (asked when loss >=10 lb).
validweight=np.where(d.cycle=='2011-2012',d.WHD020.between(70,464)&d.WHD050.between(70,485),d.WHD020.between(75,493)&d.WHD050.between(75,559)); loss=d.WHD050-d.WHD020
d['fried_loss_lb']=np.where(validweight,loss,np.nan)
d['fried_shrinking']=np.where(validweight&(loss<=10),0,np.where(validweight&(loss>10)&d.WHQ060.isin([1,2]),(d.WHQ060==2).astype(float),np.nan))
d['fried_exhaustion_proxy']=np.where(d.DPQ040.isin([0,1,2,3]),(d.DPQ040>=2).astype(float),np.nan)
d['fried_weakness_proxy']=d.fi_weak_grip;d['fried_mobility_proxy']=d.fi_walk
# CDC GPAQ MET values used to estimate energy, imperfect comparator to CHS kcal/week.
parts=[]
for q,day,minute,met in [('PAQ605','PAQ610','PAD615',8),('PAQ620','PAQ625','PAD630',4),('PAQ635','PAQ640','PAD645',4),('PAQ650','PAQ655','PAD660',8),('PAQ665','PAQ670','PAD675',4)]:
 parts.append(np.where(d[q]==2,0,np.where((d[q]==1)&d[day].between(1,7)&d[minute].between(10,1439),d[day]*d[minute]*met,np.nan)))
pa=np.array(parts).T;d['activity_met_minutes_week']=np.where(np.isfinite(pa).all(axis=1),np.sum(pa,axis=1),np.nan)
d['activity_kcal_week_proxy']=d.activity_met_minutes_week*d.BMXWT/60
d['fried_low_activity_proxy']=np.where(d.activity_kcal_week_proxy.notna()&(male|female),(d.activity_kcal_week_proxy<np.where(male,383,270)).astype(float),np.nan)
fried=['fried_shrinking','fried_exhaustion_proxy','fried_weakness_proxy','fried_mobility_proxy','fried_low_activity_proxy'];d['fried_adapted_count']=d[fried].sum(axis=1,min_count=5);d['fried_adapted_3plus']=np.where(d.fried_adapted_count.notna(),(d.fried_adapted_count>=3).astype(float),np.nan)
phys=[x for x in fried if x!='fried_exhaustion_proxy'];d['fried_physical4']=d[phys].mean(axis=1).where(d[phys].notna().all(axis=1))
# Recreate earlier four-item score with original definitions for a same-participant comparison.
oldweak=np.where((d.MGDEXSTS==1)&d.MGDCGSZ.notna(),np.where(male,d.MGDCGSZ<52,d.MGDCGSZ<32),np.nan);a=['PAQ605','PAQ620','PAQ635','PAQ650','PAQ665'];oldactivity=np.where(d[a].isin([1,2]).all(axis=1),(d[a]==2).all(axis=1).astype(float),np.nan);oldloss=np.where(validweight,(loss/d.WHD050>=.05).astype(float),np.nan)
d['prior_physical4']=pd.DataFrame({'w':oldweak,'m':d.fi_walk,'a':oldactivity,'l':oldloss}).mean(axis=1).where(pd.DataFrame({'w':oldweak,'m':d.fi_walk,'a':oldactivity,'l':oldloss}).notna().all(axis=1))
d['sleep_told_doctor']=binary_yes(d.SLQ050);d['ever_smoked']=binary_yes(d.SMQ020);d['education_clean']=numeric_valid(d.DMDEDUC2,[1,2,3,4,5]);d['wtmec4yr']=d.WTMEC2YR/2
older=d.loc[d.RIDAGEYR>=60].copy();older.to_csv(P/'nhanes_2011_2014_clean.csv',index=False)
scorecols=['cycle','SEQN','RIDAGEYR','RIAGENDR','SDMVSTRA','SDMVPSU','WTMEC2YR','wtmec4yr','sleep_told_doctor']+items+[c for n in versions for c in [n,n+'_n',n+'_complete']]+fried+['dominant_grip_mean','grip_threshold_kg','fried_loss_lb','activity_met_minutes_week','activity_kcal_week_proxy','fried_adapted_count','fried_adapted_3plus','fried_physical4','prior_physical4'];older[scorecols].to_csv(P/'frailty_scores.csv',index=False)
pd.DataFrame(itemdict).to_csv(P/'fi_item_dictionary.csv',index=False)
# Audit observed components in age-60+ MEC respondents by cycle: unweighted, descriptive.
audit=[]
for cyc,g in older.loc[older.WTMEC2YR>0].groupby('cycle'):
 for v in items:
  vals=g[v];ok=vals.notna();rho=stats.spearmanr(g.loc[ok,'RIDAGEYR'],vals[ok]).statistic if ok.sum()>2 else np.nan
  audit.append({'cycle':cyc,'item':v,'n_eligible':len(g),'n_observed':int(ok.sum()),'missing_fraction':1-ok.mean(),'prevalence_any_deficit':(vals[ok]>0).mean(),'prevalence_full_deficit':(vals[ok]==1).mean(),'spearman_age':rho})
pd.DataFrame(audit).to_csv(P/'fi_item_diagnostics.csv',index=False)
pairs=[]
for i,a in enumerate(items):
 for b in items[i+1:]:
  pair=older.loc[(older.WTMEC2YR>0),[a,b]].dropna()
  if len(pair)>2 and pair[a].std()>0 and pair[b].std()>0:
   rho=pair[a].corr(pair[b]);
   if abs(rho)>.80:pairs.append({'item_a':a,'item_b':b,'n':len(pair),'pearson_r':rho})
pd.DataFrame(pairs,columns=['item_a','item_b','n','pearson_r']).to_csv(P/'fi_correlated_items.csv',index=False)
# Sampling design includes all MEC respondents; domain selection occurs after design creation.
z=d.loc[(d.WTMEC2YR>0)&d.SDMVSTRA.notna()&d.SDMVPSU.notna()].copy().reset_index(drop=True)
strata=(z.cycle.astype(str)+':'+z.SDMVSTRA.astype(int).astype(str)).to_numpy();psu=z.SDMVPSU.to_numpy();w=z.wtmec4yr.to_numpy()
df=len(set(zip(strata,psu)))-len(set(strata)); assert all(len(np.unique(psu[strata==s]))>=2 for s in np.unique(strata))
def V(u):
 u=np.atleast_2d(u);v=np.zeros((u.shape[0],u.shape[0]))
 for s in np.unique(strata):
  sel=strata==s;ids=np.unique(psu[sel]);sums=np.stack([u[:,sel&(psu==p)].sum(axis=1) for p in ids]);v+=len(ids)/(len(ids)-1)*(sums-sums.mean(axis=0)).T@(sums-sums.mean(axis=0))
 return v
rows=[]
def add(kind,label,n,b,se=np.nan,p=np.nan,adjustment='',cycle='pooled',family=''):
 localdf=df if cycle=='pooled' else len(set(zip(strata[z.cycle==cycle],psu[z.cycle==cycle])))-len(set(strata[z.cycle==cycle]));t=stats.t.ppf(.975,localdf);rows.append(dict(kind=kind,outcome=label,cycle=cycle,adjustment=adjustment,n=int(n),estimate=float(b),se=float(se),ci_low=b-t*se,ci_high=b+t*se,p_value=p,design_df=localdf,family=family))
def fit(y,mask,adj='age_sex',label=None,cycle='pooled',family='exploratory'):
 m=np.asarray(mask)&z[y].notna().to_numpy()&z.sleep_told_doctor.notna().to_numpy()
 if adj=='full':m &= z[['education_clean','INDFMPIR','ever_smoked','RIDRETH3']].notna().all(axis=1).to_numpy()
 X=np.column_stack([np.ones(len(z)),z.sleep_told_doctor.fillna(0),(z.cycle=='2013-2014').astype(int)]) if cycle=='pooled' else np.column_stack([np.ones(len(z)),z.sleep_told_doctor.fillna(0)])
 if adj!='unadjusted':X=np.column_stack([X,z.RIDAGEYR.fillna(0)-70,(z.RIAGENDR==2).astype(int)])
 if adj=='full':X=np.column_stack([X,z.INDFMPIR.fillna(0),z.ever_smoked.fillna(0),*[(z.education_clean==v).astype(int) for v in [2,3,4,5]],*[(z.RIDRETH3==v).astype(int) for v in [2,3,4,6,7]]])
 yval=z[y].fillna(0).to_numpy();xx=X[m];ww=w[m];bread=np.linalg.inv(xx.T@(ww[:,None]*xx));beta=bread@(xx.T@(ww*yval[m]));res=yval[m]-xx@beta;u=np.zeros((X.shape[1],len(z)));u[:,m]=bread@(xx.T*(ww*res));se=np.sqrt(np.diag(V(u)));localdf=df if cycle=='pooled' else len(set(zip(strata[z.cycle==cycle],psu[z.cycle==cycle])))-len(set(strata[z.cycle==cycle]));p=2*stats.t.sf(abs(beta[1]/se[1]),localdf);add('weighted_linear',label or y,m.sum(),beta[1],se[1],p,adj,cycle,family)
def mean(y,mask,label):
 m=np.asarray(mask)&z[y].notna().to_numpy();yy=z[y].fillna(0).to_numpy();ww=w*m;est=(ww*yy).sum()/ww.sum();u=ww*(yy-est)/ww.sum();add('weighted_mean',label,m.sum(),est,np.sqrt(V(u)[0,0]))
def corr(a,b,mask,label):
 m=np.asarray(mask)&z[[a,b]].notna().all(axis=1).to_numpy();x=z.loc[m,a].to_numpy();y=z.loc[m,b].to_numpy();ww=w[m];mx=np.average(x,weights=ww);my=np.average(y,weights=ww);co=np.average((x-mx)*(y-my),weights=ww)/np.sqrt(np.average((x-mx)**2,weights=ww)*np.average((y-my)**2,weights=ww));add('weighted_correlation',label,m.sum(),co)
olderdomain=(z.RIDAGEYR>=60).to_numpy();common=olderdomain&z[['sleep_told_doctor','fi_all','fi_no_sleep','fi_no_fatigue','fi_neither']].notna().all(axis=1).to_numpy();commoncomparison=common&z[['prior_physical4','fried_physical4']].notna().all(axis=1).to_numpy()
for n in ['fi_all','fi_no_sleep','fi_no_fatigue','fi_neither']:
 for adj in ['unadjusted','age_sex','full']:fit(n,common,adj,family='primary' if adj=='age_sex' else 'exploratory')
 mean(n,common,n); z[n+'_high']=(z[n]>.21).astype(float);mean(n+'_high',common,n+'_high_above_0.21')
for n in ['fi_no_sleep','fi_no_fatigue','fi_neither']:
 z['delta']=z.fi_all-z[n];fit('delta',common,'age_sex','paired_fi_all_minus_'+n,family='primary');corr('fi_all',n,common,'fi_all_vs_'+n)
for n in ['fi_all','fi_neither','prior_physical4','fried_physical4']:
 fit(n,commoncomparison,'age_sex','comparison_'+n,family='comparison');mean(n,commoncomparison,'comparison_'+n)
corr('fi_neither','prior_physical4',commoncomparison,'fi_neither_vs_prior_physical4')
for cmp in ['prior_physical4','fried_physical4']:
 z['comparison_delta']=z.fi_neither-z[cmp];fit('comparison_delta',commoncomparison,'age_sex','paired_fi_neither_minus_'+cmp,family='comparison')
for n in ['fi_all','fi_neither']:
 fit(n,common&(z.RIDAGEYR>=65),'age_sex','age65_'+n)
 fit(n+'_complete',olderdomain,'age_sex','strict_complete_'+n)
 for cyc in ['2011-2012','2013-2014']:fit(n,common&(z.cycle==cyc),'age_sex',n,cyc)
fit('fi_no_mood',common,'age_sex','fi_no_mood');fit('fi_neither_no_care',common,'age_sex','fi_neither_no_care');fit('fi_missing5_subset',olderdomain,'age_sex','fi_missing5_subset');print('strict missingness items',len(strictitems),strictitems)
for n in ['fried_adapted_count','fried_physical4']:
 fit(n,olderdomain,'age_sex',n);mean(n,olderdomain,n)
mean('fried_adapted_3plus',olderdomain,'fried_adapted_3plus_not_diagnostic')
# Unweighted distribution and missingness checks; age trends restricted to older adults.
for n in ['fi_all','fi_neither','fried_physical4','prior_physical4']:
 for group,mask in [('all',olderdomain),('male',olderdomain&(z.RIAGENDR==1)),('female',olderdomain&(z.RIAGENDR==2))]:
  g=z.loc[mask,n].dropna();add('unweighted_quantile',n+'_'+group+'_median',len(g),g.median());add('unweighted_quantile',n+'_'+group+'_p99',len(g),g.quantile(.99));mean(n,mask,n+'_'+group)
 for ages in [(60,64),(65,69),(70,74),(75,79),(80,80)]:mean(n,olderdomain&(z.RIDAGEYR>=ages[0])&(z.RIDAGEYR<=ages[1]),n+'_age_'+str(ages[0])+'_'+str(ages[1]))
r=pd.DataFrame(rows)
for fam in ['primary','comparison','exploratory']:
 idx=r.index[(r.family==fam)&r.p_value.notna()];pv=r.loc[idx,'p_value'].to_numpy();order=np.argsort(pv);corrected=np.minimum(1,np.maximum.accumulate(pv[order]*np.arange(len(pv),0,-1)));r.loc[idx[order],'p_holm_within_family']=corrected
r.to_csv(P/'survey_results.csv',index=False)
print('design df',df,'n elder',int(olderdomain.sum()),'common FI',int(common.sum()),'common comparison',int(commoncomparison.sum()))
print(r.loc[(r.kind=='weighted_linear')&(r.adjustment=='age_sex'),['outcome','cycle','n','estimate','ci_low','ci_high','p_value','p_holm_within_family']].to_string(index=False))
