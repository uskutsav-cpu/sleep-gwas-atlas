#!/usr/bin/env python3
from pathlib import Path
import sqlite3,gzip,csv,numpy as np,pandas as pd
from scipy.stats import norm
root=Path(__file__).resolve().parents[1];work=root/'work';out=root/'outputs';db=sqlite3.connect(work/'stage6_psychad_rbm6_model.db');ws={r[0]:r[1] for r in db.execute('select rsid,weight from weights')};db.close();cov={}
with gzip.open(work/'stage6_psychad_rbm6_cov.txt.gz','rt') as f:
 next(f)
 for l in f:
  gene,a,b,v=l.split();cov[(a,b)]=cov[(b,a)]=float(v)
ids=list(ws); C=np.array([[cov[(a,b)] for b in ids] for a in ids]);eigs=np.linalg.eigvalsh(C);rows=[]
for official in pd.read_csv(out/'stage6_psychad_rbm6_twas_exploratory.csv').to_dict('records'):
 name=official['trait'];g=pd.read_csv(work/f'stage6_rbm6_{name}_gwas.tsv',sep='\t').set_index('rsid');keys=list(g.index);w=np.array([ws[k] for k in keys]);c=np.array([[cov[(a,b)] for b in keys] for a in keys]);z=g.loc[keys,'beta'].to_numpy()/g.loc[keys,'se'].to_numpy();var=float(w@c@w);manual=float(np.sum(w*np.sqrt(np.diag(c))*z)/np.sqrt(var));rows.append({'trait':name,'manual_z':manual,'official_z':official['zscore'],'difference':abs(manual-official['zscore']),'manual_p':2*norm.sf(abs(manual)),'used_variance':var,'full_covariance_min_eigenvalue':float(eigs.min())})
assert max(r['difference'] for r in rows)<1e-6
with open(out/'stage6_psychad_rbm6_formula_check.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('max abs z difference',max(r['difference'] for r in rows),'minimum covariance eigenvalue',eigs.min())
