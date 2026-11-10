#!/usr/bin/env python3
from pathlib import Path
import json,pandas as pd
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';j=json.loads((o/'stage6_foxp2_ukb_hapmap3_candidate_metadata.json').read_text());x=pd.read_csv(w/'stage6_foxp2_susie_aligned.csv').set_index('rs');allowed=set((w/'stage6_foxp2_fin_common_ids.txt').read_text().splitlines());lift=LiftOver(str(w/'hg38ToHg19.over.chain.gz')); rows=[]
for r in j['matches']:
 rs=r['SNP']
 if rs not in allowed:continue
 d=x.loc[rs];m=lift.convert_coordinate('chr7',int(d.pos)-1)
 if len(m)!=1 or m[0][0]!='chr7' or m[0][2]!='+' or m[0][1]+1!=int(r['POS']):continue
 if {r['A1'].upper(),r['A2'].upper()}!={d.anchor_effect_allele,d.anchor_other_allele}:continue
 rows.append({'rs':rs,'pos37':int(r['POS']),'A1':r['A1'].upper(),'A2':r['A2'].upper(),'anchor':d.anchor_effect_allele,'sign':1 if r['A1'].upper()==d.anchor_effect_allele else -1})
assert len({r['rs'] for r in rows})==len(rows);(w/'stage6_foxp2_ukb_eligible.json').write_text(json.dumps(rows));print('matched, coordinate and allele checked',len(rows),'candidate present',any(r['rs']=='rs12536335' for r in rows),'coloc lead present',any(r['rs']=='rs35624276' for r in rows))
