#!/usr/bin/env python3
"""Aggregate call properties of preselected GEPITEL tags; no outcomes or sample IDs."""
import csv,hashlib,json,math
from collections import Counter
from pathlib import Path
import numpy as np
root=Path('.')
source=root/'outputs/stage6_gepitel_full_reference_tags.csv'
bim=root/'work/gepitel/GEPITEL.bim'; fam=root/'work/gepitel/GEPITEL.fam'; bed=root/'work/gepitel/GEPITEL.bed'
rows=[r for r in csv.DictReader(source.open()) if r['best_r2'] and float(r['best_r2'])>=.8]
assert len(rows)==11
idx={}
with bim.open() as f:
 for i,line in enumerate(f):
  c,rs,_,pos,a1,a2=line.split()
  for r in rows:
   names=r['best_manifest_probe'].split(';')
   if c==r['chr_GRCh37'] and pos==r['best_tag_pos'] and rs in names:
    idx.setdefault(r['variant'],[]).append((i,rs,a1,a2))
n=sum(1 for _ in fam.open()); nmarkers=i+1; step=(n+3)//4
arr=np.memmap(bed,mode='r',dtype=np.uint8)
assert bytes(arr[:3])==b'\x6c\x1b\x01' and len(arr)==3+nmarkers*step
out=[]
for r in rows:
 hits=idx.get(r['variant'],[])
 x={'candidate':r['variant'],'score_lists':r['score_lists'],'reference_tag':r['best_tag_id'],'reference_r2':r['best_r2'],'matched_bim_rows':len(hits)}
 if len(hits)==1:
  i,rs,a1,a2=hits[0]; raw=arr[3+i*step:3+(i+1)*step]; codes=((raw[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n]
  counts=Counter(map(int,codes)); assert sum(counts.values()) == n
  nmiss=counts[1]; called=n-nmiss
  x.update({'bim_marker':rs,'bim_allele1':a1,'bim_allele2':a2,'n_records':n,'n_called':called,'n_missing':nmiss,'call_rate':round(called/n,6),'n_homozygous_a1':counts[0],'n_heterozygous':counts[2],'n_homozygous_a2':counts[3], 'maf_in_called':round(min((2*counts[0]+counts[2])/(2*called),(2*counts[3]+counts[2])/(2*called)),6) if called else ''})
 out.append(x)
fields=['candidate','score_lists','reference_tag','reference_r2','matched_bim_rows','bim_marker','bim_allele1','bim_allele2','n_records','n_called','n_missing','call_rate','n_homozygous_a1','n_heterozygous','n_homozygous_a2','maf_in_called']
with open('outputs/stage6_gepitel_reference_tag_call_summary.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
 return h.hexdigest()
receipt={'n_fam_records':n,'selected_candidates':len(rows),'unique_bim_matches':sum(x['matched_bim_rows']==1 for x in out),'input_sha256':{str(p):sha(p) for p in [source,bim,fam,bed]},'missing_call_range':[min(x['n_missing'] for x in out),max(x['n_missing'] for x in out)],'maf_range':[min(x['maf_in_called'] for x in out),max(x['maf_in_called'] for x in out)],'caveat':'Target array genotype call summaries only. No candidate genotype, cohort-specific proxy r2, phase, phenotype linkage, prediction, or score validation. IDs of participants were not extracted.'}
Path('outputs/stage6_gepitel_reference_tag_call_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
