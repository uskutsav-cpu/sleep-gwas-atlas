#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
from pyliftover import LiftOver
O=Path('outputs'); W=Path('work')
q=pd.read_csv(O/'stage6_rbm6_gtex_cortex_region.tsv',sep='\t');q=q[(q.type=='SNP') & (q.ref.str.len()==1) & (q.alt.str.len()==1)].copy();q['pair']=q.apply(lambda r:''.join(sorted([r.ref,r.alt])),axis=1);q=q[~q.pair.isin(['AT','CG']) & q.se.gt(0)].drop_duplicates('variant');assert not q.duplicated(['position','pair']).any()
x=q
for name in ['Mak_FI','Ye_FFS']:
 d=pd.read_csv(W/f'stage3_inka1_{name}_aligned.csv');d['pair']=d.apply(lambda r:''.join(sorted([r.ref,r.alt])),axis=1);d=d.rename(columns={'position_GRCh38':'position','ref':'gref','alt':'galt','gwas_beta_alt':name+'_beta','gwas_se':name+'_se','gwas_p':name+'_p','gwas_af_alt':name+'_af'});d=d[['position','pair','gref','galt',name+'_beta',name+'_se',name+'_p',name+'_af']];x=x.merge(d,on=['position','pair'],validate='one_to_one');assert ((x.ref==x.gref)&(x.alt==x.galt)).all();x=x.drop(columns=['gref','galt'])
lo=LiftOver(str(W/'hg38ToHg19.over.chain.gz'));defpos=lambda p:lo.convert_coordinate('chr3',int(p)-1)
positions=[]
for p in x.position:
 m=defpos(p);positions.append(m[0][1]+1 if len(m)==1 and m[0][0]=='chr3' and m[0][2]=='+' else np.nan)
x['pos37']=positions;x=x.dropna(subset=['pos37']);x['pos']=x.position;x['rs']=x.variant;x['anchor_effect_allele']=x.alt;x['anchor_other_allele']=x.ref
assert x.rs.is_unique and not x.duplicated(['pos37','pair']).any();x.to_csv(W/'stage6_rbm6_gtex_direct_input.csv',index=False);print('both-trait input',len(x),'top variant present',x.variant.eq('chr3_50085532_T_C').any(),flush=True)
