from pathlib import Path
import pandas as pd,json
from pyliftover import LiftOver
O=Path('outputs');W=Path('work');d=pd.read_csv(O/'stage6_rbm6_five_gwas_common_region.csv.gz');l=LiftOver(str(W/'hg38ToHg19.over.chain.gz'));pos=[]
for x in d.pos:
 a=l.convert_coordinate('chr3',int(x)-1);pos.append(a[0][1]+1 if len(a)==1 and a[0][0]=='chr3' and a[0][2]=='+' else None)
d['pos37']=pos;d=d[d.pos37.between(49750000,50750000)].copy();d['rs']=d.rsid;d['anchor_effect_allele']=d.alt;d['anchor_other_allele']=d.ref;assert d.rs.is_unique and not d.duplicated(['pos','pair']).any();d.to_csv(W/'stage6_rbm6_gnomad_pair_input.csv',index=False);print('selected',len(d),'GRCh37 min/max',int(d.pos37.min()),int(d.pos37.max()),'HFRS lead',bool(d.rs.eq('rs2282749').any()),'Jansen original lead',bool(d.rs.eq('rs10865954').any()))
