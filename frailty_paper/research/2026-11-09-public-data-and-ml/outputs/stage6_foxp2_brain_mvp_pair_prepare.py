from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import norm
from pyliftover import LiftOver
O=Path('outputs'); W=Path('work')
q=pd.read_csv(O/'stage6_singlebrain_ext_foxp2_full_regional.tsv.gz',sep='\t'); q=q[(q.fixed_sd>0)&q.Allele.isin(['A','C','G','T'])&(q.Allele.eq(q.ref)|q.Allele.eq(q.alt))].copy(); q['pair']=q.apply(lambda r:''.join(sorted([str(r.ref),str(r.alt)])),axis=1)
m=pd.read_csv(W/'stage6_foxp2_MVP_insomnia_region.tsv.gz',sep='\t'); m=m[m.effect_allele.isin(['A','C','G','T'])&m.other_allele.isin(['A','C','G','T'])&(m.odds_ratio>0)&m.p_value.between(0,1,inclusive='neither')].copy(); m['pair']=m.apply(lambda r:''.join(sorted([r.effect_allele,r.other_allele])),axis=1)
x=q.merge(m,left_on=['pos','pair'],right_on=['base_pair_location','pair']); x=x[~x.pair.isin(['AT','CG'])].copy();x['anchor_effect_allele']=x.effect_allele;x['anchor_other_allele']=x.other_allele;x['rs']=x.rsid;x['brain_z_anchor']=np.where(x.Allele==x.effect_allele,x.fixed_beta,-x.fixed_beta)/x.fixed_sd;x['mvp_z_anchor']=np.sign(np.log(x.odds_ratio))*norm.isf(x.p_value/2)
lift=LiftOver(str(W/'hg38ToHg19.over.chain.gz'));v=x.pos.map(lambda p:lift.convert_coordinate('chr7',int(p)-1));x['pos37']=v.map(lambda a:a[0][1]+1 if len(a)==1 and a[0][0]=='chr7' and a[0][2]=='+' else np.nan);x=x[x.pos37.between(113543159,114543159)].copy(); assert x.rs.is_unique and not x.duplicated(['pos','pair']).any();x.to_csv(W/'stage6_foxp2_brain_mvp_pair_input.csv',index=False);print('paired nonpalindromic SNPs',len(x))
