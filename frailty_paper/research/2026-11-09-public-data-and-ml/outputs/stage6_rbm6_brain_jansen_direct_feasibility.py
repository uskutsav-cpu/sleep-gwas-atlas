from pathlib import Path
import pandas as pd,json
from pyliftover import LiftOver
O=Path('outputs');W=Path('work');q=pd.read_csv(O/'stage6_singlebrain_ext_rbm6_full_regional.tsv.gz',sep='\t');q=q[(q.fixed_sd>0)&(q.Allele.eq(q.ref)|q.Allele.eq(q.alt))].copy();l=LiftOver(str(W/'hg38ToHg19.over.chain.gz'))
def conv(p):
 a=l.convert_coordinate('chr3',int(p)-1);return a[0][1]+1 if len(a)==1 and a[0][0]=='chr3' and a[0][2]=='+' else None
q['pos37']=q.pos.map(conv);q=q[q.pos37.between(49750000,50750000)].copy();q['pair']=q.apply(lambda r:''.join(sorted([str(r.ref),str(r.alt)])),axis=1);q=q[~q.pair.isin(['AT','CG'])];j=pd.read_csv(W/'stage3_region_cache/Jansen_insomnia.regional.tsv.gz',sep='\t');j=j[j.CHR.astype(str).eq('3')].copy();j['pair']=j.apply(lambda r:''.join(sorted([str(r.A1),str(r.A2)])),axis=1);x=q.merge(j,left_on=['pos37','pair'],right_on=['BP','pair']);receipt={'qtl_valid_in_window':len(q),'exact_position_allele_matches_before_reference':len(x),'unique_jansen_ids':x.SNP.nunique(),'five_trait_reference_qtl_matches':541,'limitation':'Feasibility only: no reference LD matching, colocalization or source variant-ID validation'};(O/'stage6_rbm6_brain_jansen_direct_feasibility.json').write_text(json.dumps(receipt,indent=2)+'\n');print(receipt)
