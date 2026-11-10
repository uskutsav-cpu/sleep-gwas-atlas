#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,json
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';b=np.load(w/'stage6_foxp2_1000g_eur_genotypes.npz');lookup={x:i for i,x in enumerate(b['rsid'])};fin_idx=np.loadtxt(w/'stage6_foxp2_fin_reference_indices.txt',dtype=int)-1;reports=[]
for name in ['left','right']:
 d=pd.read_csv(o/f'stage6_foxp2_gnomad_fin_{name}_variants.csv');n=len(d);R=np.fromfile(w/f'stage6_foxp2_gnomad_fin_{name}_ld.bin',dtype='<f8').reshape(n,n);g=b['dosage'][[lookup[x] for x in d.rs]].astype(np.float64);g-=g.mean(axis=1,keepdims=True);g/=np.linalg.norm(g,axis=1,keepdims=True);eur=g@g.T;ii,jj=np.where(np.triu((abs(R)>.1)&(abs(eur)>.1),1));same=(R[ii,jj]*eur[ii,jj])>0
 k=list(d.rs);a=k.index('rs35624276');c=k.index('rs12536335');print(name,'anchor signed r',R[a,c],'r2',R[a,c]**2,'pairs',len(ii),'opposite',int((~same).sum()),'corr',np.corrcoef(R[ii,jj],eur[ii,jj])[0,1],flush=True);reports.append({'window':name,'n_variants':n,'anchor_r':float(R[a,c]),'anchor_r2':float(R[a,c]**2),'n_pairs_both_abs_r_gt_01':len(ii),'opposite_sign_pairs':int((~same).sum()),'signed_r_correlation':float(np.corrcoef(R[ii,jj],eur[ii,jj])[0,1])})
(o/'stage6_foxp2_gnomad_fin_reference_check.json').write_text(json.dumps(reports,indent=2)+'\n')
