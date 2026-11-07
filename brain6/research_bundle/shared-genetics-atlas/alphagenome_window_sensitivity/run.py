#!/usr/bin/env python3
"""Descriptive re-summaries of saved, user-authorized AlphaGenome tracks; no model training."""
import csv,gzip,hashlib,json,math,struct,collections
from pathlib import Path
root=Path('outputs/brain6-alphagenome-pilot/results/jbrowse_expanded'); out=Path('outputs/shared-genetics-atlas/alphagenome_window_sensitivity')
with (root/'expanded_variant_evidence.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
results=[]
for row in rows:
 assert row['tissue_prediction_status']=='PREDICTED'
 folder=root/'variants'/hashlib.sha256(row['atlas_variant'].encode()).hexdigest()[:20]
 m=json.loads((folder/'provenance.json').read_text()); raw=(folder/'track_window.f32.gz').read_bytes()
 assert hashlib.sha256(raw).hexdigest()==m['track_window_sha256']==row['prediction_window_sha256']
 a=struct.unpack('<16008f',gzip.decompress(raw))
 tracks=[i for i,t in enumerate(m['tracks']) if t['output_type']=='dnase' and t['metadata'].get('biosample_life_stage')=='adult'];assert len(tracks)==1
 i=tracks[0];ref=a[i*4002:i*4002+2001];alt=a[i*4002+2001:(i+1)*4002]
 r={'region':row['candidate_locus_id'],'snp':row['snp'],'variant':row['atlas_variant'],'input_sha256':m['track_window_sha256']}
 for width in [101,501,2001]:
  lo=1000-width//2; hi=lo+width; rs=sum(ref[lo:hi]); als=sum(alt[lo:hi]);assert rs>=0 and als>=0
  r[f'ref_{width}']=rs;r[f'alt_{width}']=als
  for pseudo in [0.1,1,10]:r[f'log2ratio_{width}_pc{pseudo}']=math.log2((als+pseudo)/(rs+pseudo))
 assert abs(r['log2ratio_501_pc1']-float(row['dnase_log2_ratio_501bp']))<1e-10
 results.append(r)
with (out/'window_scores.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
print('variants',len(results))
for width in [101,501,2001]:
 for pseudo in [0.1,1,10]:
  k=f'log2ratio_{width}_pc{pseudo}'; ordered=sorted(results,key=lambda x:abs(x[k]),reverse=True);idx=next(i for i,r in enumerate(ordered) if r['snp']=='rs62264764');print(width,pseudo,'top',ordered[0]['snp'],ordered[0][k],'rs62264764 rank',idx+1,'score',ordered[idx][k])
