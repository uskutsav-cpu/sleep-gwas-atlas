#!/usr/bin/env python3
"""Join previously selected sleep hit models to public ADHD model outputs."""
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;parent=p.parent
d=json.loads((p/'receipt.json').read_text());bykey={}
for x in d['matches']:
 name=Path(x['member']).name.split('.');assert name[0]=='fusion'
 tissue,mode=name[1:3];cols=x['header'].split('\t');vals=x['row'].split('\t');assert len(cols)==len(vals)
 row=dict(zip(cols,vals));key=(tissue,mode,row['ID']);bykey.setdefault(key,[]).append(row)
assert len(d['regular_members'])==400 and len(d['matches'])==165
out=[];same_model=['CHR','P0','P1','HSQ','EQTL.ID','EQTL.R2','EQTL.Z','NSNP','NWGT','MODEL','MODELCV.R2','MODELCV.PV']
for source,mode in (('twas_hits.gtex-ddp_selected.tsv','latent_full'),('twas_hits.gtex-rddp_selected.tsv','latent_residual')):
 with (parent/source).open() as f:
  for s in csv.DictReader(f,delimiter='\t'):
   if s['TRAIT']!='UKB_1200_Sleeplessness_or_insomnia':continue
   k=(s['TISSUE'],mode,s['ID']);matches=bykey.get(k,[]);assert len(matches)==1
   a=matches[0];assert all(s[x]==a[x] for x in same_model)
   sp,ap=float(s['TWAS.P']),float(a['TWAS.P']);assert math.isfinite(sp) and math.isfinite(ap)
   out.append({'tissue':s['TISSUE'],'model':mode,'feature':s['ID'],'sleep_z':s['TWAS.Z'],'sleep_p':s['TWAS.P'],
    'sleep_coloc_PP4':s['COLOC.PP4'],'adhd_z':a['TWAS.Z'],'adhd_p':a['TWAS.P'],
    'adhd_coloc_PP4':a['COLOC.PP4'],'model':mode,'weights_model':s['MODEL'],'weight_count':s['NWGT'],
    'molecular_lead':s['EQTL.ID'],'sleep_gwas_lead':s['BEST.GWAS.ID'],'adhd_gwas_lead':a['BEST.GWAS.ID']})
assert len(out)==8 and sum(float(x['adhd_p'])<.05 for x in out)==1
target=p/'matched_models.tsv'
with target.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys(),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'archive_sha256':d['sha256'],'archive_bytes':d['bytes_read'],'models_matched':len(out),'source_model_fields_checked':same_model,
 'adhd_nominal_p_under_0_05':sum(float(x['adhd_p'])<.05 for x in out),
 'smallest_adhd_p':min(float(x['adhd_p']) for x in out),'adhd_pp4_all_unavailable':all(x['adhd_coloc_PP4']=='NA' for x in out),
 'matched_table_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
 'note':'Retrospective selected-model comparison; matching model descriptors is not proof of independent replication or shared causal variants.'}
(p/'comparison_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
