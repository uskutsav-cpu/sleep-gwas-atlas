#!/usr/bin/env python3
import csv,gzip,json,hashlib
from pathlib import Path
ids=('rs2269457','rs2585408','rs10462028')
ref={}
for rs in ids:
 d=json.load(open('work/hunt_'+rs+'_ensembl38.json'));maps=[x for x in d['mappings'] if x.get('assembly_name')=='GRCh38' and x.get('seq_region_name') in [str(i) for i in range(1,23)]]
 if len(maps)!=1: raise ValueError('not uniquely mapped '+rs)
 ref[rs]=(str(maps[0]['seq_region_name']),str(maps[0]['start']),set(maps[0]['allele_string'].split('/')))
p=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/finngen_R13_F5_INSOMNIA.gz'); results=[];n=0
with gzip.open(p,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  n+=1
  for rs,(ch,pos,alleles) in ref.items():
   if row['#chrom']==ch and row['pos']==pos:
    results.append(dict(rsid=rs,chr38=ch,pos38=pos,source_id=row['rsids'],ref=row['ref'],alt=row['alt'],reference_pair_match=set([row['ref'],row['alt']])==alleles,beta_alt=row['beta'],se=row['sebeta'],p=row['pval'],af_alt=row['af_alt']))
h=hashlib.sha256()
with open(p,'rb') as f:
 for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
with open('outputs/stage6_hunt3_finngen_insomnia_screen.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['rsid','chr38','pos38','source_id','ref','alt','reference_pair_match','beta_alt','se','p','af_alt']);w.writeheader();w.writerows(results)
print(json.dumps({'rows':n,'sha256':h.hexdigest(),'hits':results},indent=2))
