#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, json
out=Path(__file__).resolve().parent
x=pd.read_csv(out/'stage5_bryois_exc_smad5_associations.csv'); p=pd.read_csv(out/'stage5_bryois_exc_smad5_variant_positions.csv');assert x.snp_id.is_unique and p.snp_id.is_unique
x=x.merge(p,on='snp_id',validate='one_to_one');x['chr']=x.snp_id_hg38.str.split(':').str[0];x['pos']=x.snp_id_hg38.str.split(':').str[1].astype(int)
q=pd.read_csv(out/'stage5_score_queue.csv');q=q[(q.chr=='chr5')];m=x.merge(q,left_on=['chr','pos'],right_on=['chr','position_GRCh38'],how='inner',validate='many_to_many');m=m[(m.effect_allele==m.alternate_allele)&(m.other_allele==m.reference_allele)|(m.effect_allele==m.reference_allele)&(m.other_allele==m.alternate_allele)].copy();m['source_beta_alt']=np.where(m.effect_allele==m.alternate_allele,m.beta,-m.beta);m['palindromic']=m.apply(lambda r: {r.reference_allele,r.alternate_allele} in ({'A','T'},{'C','G'}),axis=1)
g=pd.read_csv(out/'stage5_variant_multitrait_gwas_matrix.csv');m=m.merge(g,left_on='variant_key',right_on='variant_GRCh38_ref_alt',how='left',validate='many_to_one');m.to_csv(out/'stage5_bryois_exc_smad5_queue_matches.csv',index=False)
print('exact position and allele-pair matches',len(m),'non-palindromic',int((~m.palindromic).sum()))
print(m[m.variant_key.isin(['chr5:136202653:A>C','chr5:136207208:G>A','chr5:136159634:A>G'])][['variant_key','snp_id','nominal_p','beta','effect_allele','source_beta_alt','Jansen_insomnia_p','Mak_FI_p']].to_string(index=False))
print('source associations in SingleBrain source-file gap',int(x.pos.between(136164133,139164640).sum()))
