#!/usr/bin/env python3
"""Get signed LD for one GRCh37 indel against existing SNPs from public gnomAD blocks."""
from pathlib import Path
import sys,requests,numpy as np,pandas as pd,json,hashlib
root=Path(__file__).resolve().parents[1];w=root/'work';o=root/'outputs';sys.path.insert(0,str(w/'ldcov'))
from ldcov.io.blockmatrix.partitioner import GridPartitioner
from ldcov.io.blockmatrix.block_codec import decode_block
fin=pd.read_csv(o/'stage6_foxp2_gnomad_narrow_variants.csv');ids=fin.rs.tolist();assert len(ids)==1257;receipts=[]
for pop in ['fin','nfe']:
 d=pd.read_csv(o/f'stage6_foxp2_gnomad_{pop}_variant_index.csv').set_index('rs').loc[ids];ind=pd.read_csv(o/f'stage6_foxp2_gnomad_{pop}_indel_index.csv').iloc[0];assert ind.ref37=='T' and ind.alt37=='TA' and ind.anchor_effect_allele=='T'
 base=f'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.adj.ld.bm';r=requests.get(base+'/metadata.json',timeout=45);r.raise_for_status();meta=r.json();bs=meta['blockSize'];part=GridPartitioner(bs,meta['nRows'],meta['nCols'],meta.get('maybeFiltered'));ix=d.gnomad_index.to_numpy();sgn=d.sign_to_anchor.to_numpy();anchor=int(ind.gnomad_index);ab=anchor//bs;vec=np.full(len(d),np.nan)
 for b in sorted(set(ix//bs)):
  lo,hi=sorted([int(b),int(ab)]);slot=part.part_slot(lo,hi)
  if slot is None:raise ValueError('required block absent, not assuming LD=0')
  name=meta['partFiles'][slot];r=requests.get(base+'/parts/'+name,timeout=180);r.raise_for_status();raw=r.content;B=decode_block(raw);sel=np.flatnonzero(ix//bs==b);loc=ix[sel]%bs;anchorloc=anchor%bs
  if b<ab:vals=B[loc,anchorloc]
  elif b>ab:vals=B[anchorloc,loc]
  else:vals=B[np.minimum(loc,anchorloc),np.maximum(loc,anchorloc)]
  vec[sel]=vals*sgn[sel]*int(ind.sign_to_anchor);receipts.append({'population':pop,'block':[lo,hi],'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'file':name});print(pop,lo,hi,'bytes',len(raw),flush=True);del raw,r,B
 assert np.isfinite(vec).all() and np.max(abs(vec))<=1.000001;pd.DataFrame({'rs':ids,'indel_anchor_r':vec}).to_csv(o/f'stage6_foxp2_gnomad_{pop}_indel_ld.csv',index=False)
(o/'stage6_foxp2_gnomad_indel_ld_receipt.json').write_text(json.dumps({'variant_GRCh38':'chr7:114534306:T>TA','effect_anchor':'T','source_blocks':receipts,'whole_matrices_checksum_verified':False},indent=2)+'\n')
