#!/usr/bin/env python3
import csv,gzip,json,hashlib,math
from pathlib import Path
src=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz');O=Path('outputs');ref=list(csv.DictReader(open(O/'stage6_foxp2_gnomad_narrow_variants.csv')));byid={r['rs']:r for r in ref};hits={k:[] for k in byid};n=0
with gzip.open(src,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  n+=1;r=byid.get(x['SNP'])
  if r and x['CHR']=='7' and x['BP']==r['pos37'] and {x['EA'],x['OA']}=={r['anchor_effect_allele'],r['anchor_other_allele']}:
   try:
    sign=1 if x['EA']==r['anchor_effect_allele'] else -1;ukb=sign*float(x['BETA_UKB1']);seu=float(x['SE_UKB1']);tg=-sign*float(x['BETA_TG']);setg=float(x['SE_TG']); assert seu>0 and setg>0 and all(map(math.isfinite,[ukb,seu,tg,setg]));hits[x['SNP']].append({'rs':x['SNP'],'pos37':r['pos37'],'anchor_effect_allele':r['anchor_effect_allele'],'source_EA':x['EA'],'source_OA':x['OA'],'ukb_beta_anchor':ukb,'ukb_se':seu,'ukb_z':ukb/seu,'tg_beta_raw':float(x['BETA_TG']),'tg_beta_anchor_adapter':tg,'tg_se':setg,'tg_z_adapter':tg/setg,'ukb_p':x['P_UKB1'],'tg_p':x['P_TG']})
   except (ValueError,AssertionError,KeyError):pass
  if n%5000000==0:print('scanned',n,flush=True)
rows=[v[0] for v in hits.values() if len(v)==1];duplicate=[k for k,v in hits.items() if len(v)>1];fields=['rs','pos37','anchor_effect_allele','source_EA','source_OA','ukb_beta_anchor','ukb_se','ukb_z','tg_beta_raw','tg_beta_anchor_adapter','tg_se','tg_z_adapter','ukb_p','tg_p'];out=O/'stage6_foxp2_fi_cohort_region_aligned.csv'
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
rec={'source':str(src),'source_was_public_archive':True,'source_whole_checksum_verified_in_this_pass':False,'rows_scanned':n,'reference_variants':len(ref),'matched_valid_unique':len(rows),'ambiguous_duplicates_excluded':len(duplicate),'twingene_sign_adapter':'negates raw source beta before anchor alignment; inferred, not author-confirmed','sha256_output':hashlib.sha256(out.read_bytes()).hexdigest()};(O/'stage6_foxp2_fi_cohort_region_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2),flush=True)
