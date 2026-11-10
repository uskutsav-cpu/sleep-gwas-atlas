#!/usr/bin/env python3
"""Compare full and direct-manifest-marker partial score sums in reference data."""
import csv, json, math, hashlib
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parent.parent; out=root/'outputs'; base=root/'work/stage4_magma_reference/g1000_eur'
def records(name):
 with (out/name).open(newline='') as f:return list(csv.DictReader(f))
direct={r['variant'] for r in records('stage6_act_array_coordinate_matches.csv') if r['name_exact']=='True'}
candidate_files=['unvalidated_15_variant_score_prototype.csv','unvalidated_22_variant_score_prototype.csv','unvalidated_22_shared_direction_score.csv']
wanted=set()
for file in candidate_files:
 for r in records(file):
  ch=r.get('chr',r.get('chr_GRCh37_assumed'));pos=r.get('position_GRCh37_inferred',r.get('position_GRCh37_assumed'));wanted.add((ch,int(pos),frozenset((r['counted_allele'],r['other_allele']))))
count=0;keys={}
with open(str(base)+'.bim') as f:
 for n,line in enumerate(f):
  ch,rs,_,bp,a,b=line.split();k=(ch,int(bp),frozenset((a,b)));count=n+1
  if k in wanted: keys.setdefault(k,[]).append((n,rs,a,b))
n=sum(1 for _ in open(str(base)+'.fam'));stride=(n+3)//4;bed=np.memmap(str(base)+'.bed',dtype='uint8',mode='r');assert bytes(bed[:3])==b'\x6c\x1b\x01' and len(bed)==3+count*stride
lut=np.array([2,np.nan,1,0],dtype=float)
def dosage(i):
 b=bed[3+i*stride:3+(i+1)*stride];return lut[((b[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n]]
sets=[('15_weighted','unvalidated_15_variant_score_prototype.csv','chr','position_GRCh37_inferred',['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']),('22_weighted','unvalidated_22_variant_score_prototype.csv','chr','position_GRCh37_assumed',['FI_beta_weight','FFS_beta_weight_inferred','UKB_insomnia_beta_weight']),('22_count','unvalidated_22_shared_direction_score.csv','chr_GRCh37_assumed','position_GRCh37_assumed',['weight'])]
res=[]
for label,file,ch,pos,weights in sets:
 rows=records(file); vals=[]
 for row in rows:
  matches=keys.get((row[ch],int(row[pos]),frozenset((row['counted_allele'],row['other_allele']))),[])
  if len(matches)!=1:raise ValueError('reference ambiguity '+row['variant'])
  idx,rs,a,b=matches[0];v=dosage(idx)
  if a!=row['counted_allele']:v=2-v
  vals.append(v)
 X=np.stack(vals,axis=1);mask=np.array([r['variant'] in direct for r in rows]);valid=np.isfinite(X).all(axis=1); X=X[valid];
 if X.shape[0]<400:raise ValueError('insufficient complete reference calls')
 for wname in weights:
  w=np.array([float(r[wname]) for r in rows]);full=X@w;part=X[:,mask]@w[mask]
  if np.std(part)==0 or np.std(full)==0:raise ValueError('constant score')
  corr=float(np.corrcoef(full,part)[0,1]);res.append({'score':label,'weight':wname,'n_total_markers':len(rows),'n_manifest_matched':int(mask.sum()),'n_reference_complete':len(X),'partial_vs_full_correlation':corr,'correlation_squared':corr*corr})
with (out/'stage6_act_direct_partial_reference_scores.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(res[0]));writer.writeheader();writer.writerows(res)
receipt={'reference':'MAGMA 1000 Genomes EUR 503 reference samples; pooled sample correlation','manifest_receipt_sha256':hashlib.sha256((out/'stage6_act_array_manifest_receipt.json').read_bytes()).hexdigest(),'method':'scores on identical reference genotypes; partial omits nonmanifest markers and uses original counted alleles and weights; full vs partial Pearson correlation, not outcome prediction','limitations':'Manifest matches are position and name only, ACT manifest version, strand, post-QC availability and imputed data unknown. All direct variants assumed available on correct contrast. Excluded candidate proxies; 1000 Genomes sample differs from ACT, FFS weight coding inferred; postselected score. Correlation is not phenotype variance explained, risk, or validation.' ,'results':res}
(out/'stage6_act_direct_partial_reference_scores.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(res,indent=2))
