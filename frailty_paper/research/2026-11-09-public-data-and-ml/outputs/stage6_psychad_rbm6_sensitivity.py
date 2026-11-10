#!/usr/bin/env python3
from pathlib import Path
import sqlite3,gzip,csv,numpy as np,pandas as pd
from scipy.stats import norm
root=Path(__file__).resolve().parents[1];wdir=root/'work';out=root/'outputs';db=sqlite3.connect(wdir/'stage6_psychad_rbm6_model.db');weights={x[0]:x[1] for x in db.execute('select rsid,weight from weights')};db.close();cov={}
with gzip.open(wdir/'stage6_psychad_rbm6_cov.txt.gz','rt') as f:
 next(f)
 for line in f:
  gene,a,b,v=line.split();cov[(a,b)]=cov[(b,a)]=float(v)
results=pd.read_csv(out/'stage6_psychad_rbm6_twas_exploratory.csv'); gwas={r.trait:pd.read_csv(wdir/f'stage6_rbm6_{r.trait}_gwas.tsv',sep='\t').set_index('rsid') for r in results.itertuples()};common=set.intersection(*[set(df.index) for df in gwas.values()]);print('common SNPs',len(common));rows=[]
def calc(g,ids):
 ids=sorted(ids);wv=np.array([weights[k] for k in ids]);c=np.array([[cov[(a,b)] for b in ids] for a in ids]);z=(g.loc[ids,'beta']/g.loc[ids,'se']).to_numpy();den=float(wv@c@wv)
 if den<=0:raise ValueError('nonpositive predicted variance')
 x=float(np.sum(wv*np.sqrt(np.diag(c))*z)/np.sqrt(den));return x,2*norm.sf(abs(x)),den
for r in results.itertuples():
 df=gwas[r.trait];full=set(df.index);z,p,v=calc(df,full);assert abs(z-r.zscore)<1e-6
 cz,cp,cv=calc(df,common);rows.append({'trait':r.trait,'scenario':'seven_GWAS_common_set','omitted_snp':'','n_snps':len(common),'zscore':cz,'pvalue':cp,'predicted_variance':cv,'full_z':z})
 for snp in sorted(full):
  zz,pp,vv=calc(df,full-{snp});rows.append({'trait':r.trait,'scenario':'leave_one_out','omitted_snp':snp,'n_snps':len(full)-1,'zscore':zz,'pvalue':pp,'predicted_variance':vv,'full_z':z})
t=pd.DataFrame(rows);t.to_csv(out/'stage6_psychad_rbm6_sensitivity.csv',index=False)
for trait,x in t.groupby('trait'):
 v=x[x.scenario=='leave_one_out'];c=x[x.scenario=='seven_GWAS_common_set'].iloc[0];print(trait,'common Z',round(c.zscore,4),'LOO range',round(v.zscore.min(),4),round(v.zscore.max(),4),'sign changes',sum(v.zscore*v.full_z<0))
