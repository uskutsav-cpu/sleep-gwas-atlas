#!/usr/bin/env python3
"""Position/allele-based descriptive match to low-grip GWAS."""
from pathlib import Path
import pandas as pd,math,json
O=Path('outputs');raw=pd.read_csv(O/'stage6_grip_gwas_raw_candidates.csv');raw.columns=raw.columns.str.strip();res=[]
for score,file,chrom,pos in [('prototype15','unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred'),('shared22','unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed')]:
 d=pd.read_csv(O/file)
 for r in d.to_dict('records'):
  group=raw[raw.variant_id.eq(r['variant'])];candidate=group[(group.chromosome.astype(str)==str(r[chrom])) & (group.base_pair_location==int(r[pos]))];candidate=candidate[candidate.apply(lambda x:{str(x.effect_allele),str(x.other_allele)}=={str(r['counted_allele']),str(r['other_allele'])},axis=1)];row={'score':score,'variant':r['variant'],'position_GRCh37':int(r[pos]),'counted_allele':r['counted_allele'],'other_allele':r['other_allele'],'source_rows_for_id':len(group),'exact_rows':len(candidate)}
  if len(candidate)==1:
   x=candidate.iloc[0];sign=1 if x.effect_allele==r['counted_allele'] else -1;row.update(source_effect_allele=x.effect_allele,source_other_allele=x.other_allele,beta_counted=sign*float(x.beta),se=float(x.standard_error),p=float(x.p_value),counted_allele_frequency=float(x.effect_allele_frequency) if sign==1 else 1-float(x.effect_allele_frequency),heterogeneity_I2=x.HetISq,heterogeneity_p=x.HetPVal)
  res.append(row)
a=pd.DataFrame(res);a.to_csv(O/'stage6_grip_gwas_candidate_aligned.csv',index=False);out={}
for group,x in a.groupby('score'):
 y=x[x.exact_rows==1];out[group]={'candidates':len(x),'exact_matches':len(y),'nominal_p_lt_005':int((y.p<.05).sum()),'missing_or_ambiguous':int((x.exact_rows!=1).sum())}
(O/'stage6_grip_gwas_candidate_align_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
