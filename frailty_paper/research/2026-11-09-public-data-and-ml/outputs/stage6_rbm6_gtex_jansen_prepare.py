#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np,json
from pyliftover import LiftOver
O=Path('outputs');W=Path('work');q=pd.read_csv(O/'stage6_rbm6_gtex_cortex_region.tsv',sep='\t');q=q[(q.type=='SNP')&(q.ref.str.len()==1)&(q.alt.str.len()==1)&(q.se>0)].copy();q['pair']=q.apply(lambda r:''.join(sorted([r.ref,r.alt])),axis=1);q=q[~q.pair.isin(['AT','CG'])].copy()
for k,g in q.groupby('variant'):
 for col in ['position','ref','alt','beta','se','ac','an']:assert g[col].nunique(dropna=False)==1,(k,col)
q=q.drop_duplicates('variant');assert not q.duplicated(['position','pair']).any();lo=LiftOver(str(W/'hg38ToHg19.over.chain.gz'));pos=[]
for p in q.position:
 m=lo.convert_coordinate('chr3',int(p)-1);pos.append(m[0][1]+1 if len(m)==1 and m[0][0]=='chr3' and m[0][2]=='+' else np.nan)
q['pos37']=pos;q=q.dropna(subset=['pos37']);j=pd.read_csv(W/'stage3_region_cache/Jansen_insomnia.regional.tsv.gz',sep='\t');j=j[j.CHR.astype(str).eq('3')].copy();j['pair']=j.apply(lambda r:''.join(sorted([str(r.A1),str(r.A2)])),axis=1);j=j.rename(columns={'BP':'pos37'});j=j.drop_duplicates(['pos37','pair'],keep=False);x=q.merge(j,on=['pos37','pair'],validate='one_to_one');x=x[(x.SE>0)&x.BETA.notna()&x.FRQ.between(0,1)].copy();x['Jansen_beta']=np.where(x.A1==x.alt,x.BETA,-x.BETA);x['Jansen_af']=np.where(x.A1==x.alt,x.FRQ,1-x.FRQ);x['pos']=x.position;x['rs']=x.variant;x['anchor_effect_allele']=x.alt;x['anchor_other_allele']=x.ref;assert x.rs.is_unique and not x.duplicated(['pos37','pair']).any();x.to_csv(W/'stage6_rbm6_gtex_jansen_input.csv',index=False);print('direct input',len(x),'top',x.variant.eq('chr3_50085532_T_C').any(),flush=True)
