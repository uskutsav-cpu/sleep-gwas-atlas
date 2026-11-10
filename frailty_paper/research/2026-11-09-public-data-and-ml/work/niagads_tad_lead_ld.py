#!/usr/bin/env python3
"""Describe pairwise reference LD between three postselected regional variants."""
import csv, hashlib, json
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent.parent; ref=root/'work/stage4_magma_reference'
source=root/'work/niagads_smad5_tad_insomnia_ld_inputs.csv'
ids=['chr5:136202653:A:C','chr5:136206558:C:T','chr5:137615744:G:A']
with source.open() as f:rows={r['variant_GRCh38']:r for r in csv.DictReader(f) if r['variant_GRCh38'] in ids}
assert set(rows)==set(ids)
indices={int(rows[x]['reference_bim_index']) for x in ids};bim={}
with (ref/'g1000_eur.bim').open() as f:
 for i,line in enumerate(f):
  if i in indices:bim[i]=line.split()
assert len(bim)==3
n=503;size=(n+3)//4; bed=np.memmap(ref/'g1000_eur.bed',dtype='uint8',mode='r');assert bytes(bed[:3])==bytes([108,27,1])
g=[]
for var in ids:
 row=rows[var];index=int(row['reference_bim_index']);refallele,alt=var.split(':')[-2:];entry=bim[index]
 assert set(entry[4:6])=={refallele,alt}
 raw=bed[3+index*size:3+(index+1)*size];codes=((raw[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n]
 dosage=np.array([2.,np.nan,1.,0.])[codes];dosage=dosage if entry[4]==alt else 2-dosage
 assert not np.isnan(dosage).all();g.append(dosage)
result=[]
for i in range(3):
 for j in range(i+1,3):
  both=np.isfinite(g[i])&np.isfinite(g[j]); r=float(np.corrcoef(g[i][both],g[j][both])[0,1])
  result.append({'variant_1':ids[i],'variant_2':ids[j],'complete_pairs':int(both.sum()),'r_alt_dosage':r,'r_squared':r*r})
receipt={'reference':'1000 Genomes phase 3 EUR, 503 samples','input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'pairs':result,'limitations':'Postselected reference genotype correlations; not Finnish or ROSMAP LD, not target causality or score validation.'}
(root/'outputs/stage6_niagads_tad_lead_reference_ld.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(result,indent=2))
