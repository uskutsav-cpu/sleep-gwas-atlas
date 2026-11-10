#!/usr/bin/env python3
"""Extract two complete windows from the public Finnish gnomAD Hail block matrix."""
from pathlib import Path
import sys,json,hashlib,requests,numpy as np,pandas as pd
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';sys.path.insert(0,str(w/'ldcov'))
from ldcov.io.blockmatrix.partitioner import GridPartitioner
from ldcov.io.blockmatrix.block_codec import decode_block
base='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.fin.common.adj.ld.bm';mr=requests.get(base+'/metadata.json',timeout=45);mr.raise_for_status();meta=mr.json();bs=meta['blockSize'];assert bs==4096;p=GridPartitioner(bs,meta['nRows'],meta['nCols'],meta.get('maybeFiltered'));choice=sys.argv[1] if len(sys.argv)>1 else 'right';assert choice in ['left','right'];selected_blocks=[1594,1595] if choice=='left' else [1595,1596];d=pd.read_csv(o/'stage6_foxp2_gnomad_fin_variant_index.csv');d=d[(d.gnomad_index//bs).isin(selected_blocks)].sort_values('gnomad_index').reset_index(drop=True);N=len(d);indices=d.gnomad_index.to_numpy();block_ids=indices//bs;groups=sorted(set(block_ids));assert groups==selected_blocks;M=np.full((N,N),np.nan,dtype=np.float64);receipt=[]
for i,a in enumerate(groups):
 for b in groups[i:]:
  slot=p.part_slot(a,b)
  if slot is None:receipt.append({'block':[int(a),int(b)],'status':'not_stored'});continue
  name=meta['partFiles'][slot];r=requests.get(base+'/parts/'+name,timeout=180);r.raise_for_status();raw=r.content;h=hashlib.sha256(raw).hexdigest();B=decode_block(raw);ii=np.flatnonzero(block_ids==a);jj=np.flatnonzero(block_ids==b);sub=B[np.ix_(indices[ii]%bs,indices[jj]%bs)].copy()
  if a==b:
   sub=sub+sub.T;np.fill_diagonal(sub,np.diag(B)[indices[ii]%bs])
  M[np.ix_(ii,jj)]=sub
  if a!=b:M[np.ix_(jj,ii)]=sub.T
  receipt.append({'block':[int(a),int(b)],'status':'retrieved','name':name,'bytes':len(raw),'sha256':h});print('read',a,b,'bytes',len(raw),flush=True);del raw,r,B,sub
sign=d.sign_to_anchor.to_numpy();M*=sign[:,None]*sign[None,:]
for name,blocks in [(choice,selected_blocks)]:
 sel=np.flatnonzero(np.isin(block_ids,blocks));x=M[np.ix_(sel,sel)];labels=d.iloc[sel].copy()
 if not np.isfinite(x).all():raise ValueError('unexpected missing correlations in '+name)
 if not np.allclose(x,x.T,atol=1e-7) or not np.allclose(np.diag(x),1,atol=1e-7):raise ValueError('not a correlation matrix')
 x.astype('<f8').tofile(w/f'stage6_foxp2_gnomad_fin_{name}_ld.bin');labels.to_csv(o/f'stage6_foxp2_gnomad_fin_{name}_variants.csv',index=False);print(name,'variants',len(sel),'range37',int(labels.pos37.min()),int(labels.pos37.max()),flush=True)
info={'matrix_url':base,'metadata_sha256':hashlib.sha256(mr.content).hexdigest(),'whole_matrix_checksum_verified':False,'blocks':receipt,'n_input':N,'unstored_block_1594_1596_not_filled':True,'ldcov_commit':'recorded separately'};(o/f'stage6_foxp2_gnomad_fin_{choice}_ld_receipt.json').write_text(json.dumps(info,indent=2)+'\n')
