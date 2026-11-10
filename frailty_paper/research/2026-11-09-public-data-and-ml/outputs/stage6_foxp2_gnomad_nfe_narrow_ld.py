#!/usr/bin/env python3
from pathlib import Path
import sys,json,hashlib,requests,numpy as np,pandas as pd
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';sys.path.insert(0,str(w/'ldcov'))
from ldcov.io.blockmatrix.partitioner import GridPartitioner
from ldcov.io.blockmatrix.block_codec import decode_block
fin=pd.read_csv(o/'stage6_foxp2_gnomad_fin_left_variants.csv');sel=fin.pos37.between(113543159,114543159);fin=fin[sel].reset_index(drop=True);n=len(fin);assert n==1257
alln=pd.read_csv(o/'stage6_foxp2_gnomad_nfe_variant_index.csv').set_index('rs');d=alln.loc[fin.rs].reset_index();assert len(d)==n
base='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm';mr=requests.get(base+'/metadata.json',timeout=40);mr.raise_for_status();meta=mr.json();bs=meta['blockSize'];p=GridPartitioner(bs,meta['nRows'],meta['nCols'],meta.get('maybeFiltered'));idx=d.gnomad_index.to_numpy();bid=idx//bs;groups=sorted(set(bid));M=np.full((n,n),np.nan);receipt=[]
for ii,a in enumerate(groups):
 for b in groups[ii:]:
  slot=p.part_slot(int(a),int(b));
  if slot is None:raise ValueError('missing block cannot be assumed zero')
  filename=meta['partFiles'][slot];r=requests.get(base+'/parts/'+filename,timeout=180);r.raise_for_status();raw=r.content;B=decode_block(raw);i=np.flatnonzero(bid==a);j=np.flatnonzero(bid==b);sub=B[np.ix_(idx[i]%bs,idx[j]%bs)].copy()
  if a==b:sub=sub+sub.T;np.fill_diagonal(sub,np.diag(B)[idx[i]%bs])
  M[np.ix_(i,j)]=sub
  if a!=b:M[np.ix_(j,i)]=sub.T
  receipt.append({'block':[int(a),int(b)],'name':filename,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()});print('read',a,b,'bytes',len(raw),flush=True);del B,raw,r,sub
sgn=d.sign_to_anchor.to_numpy();M*=sgn[:,None]*sgn[None,:];assert np.isfinite(M).all() and np.allclose(M,M.T,atol=1e-7) and np.allclose(np.diag(M),1,atol=1e-7);M.astype('<f8').tofile(w/'stage6_foxp2_gnomad_nfe_narrow_ld.bin')
# Slice the existing Finnish matrix on exactly the same ordered variants.
full=pd.read_csv(o/'stage6_foxp2_gnomad_fin_left_variants.csv');p0=len(full);F=np.fromfile(w/'stage6_foxp2_gnomad_fin_left_ld.bin',dtype='<f8').reshape(p0,p0);idxs=np.flatnonzero(sel.to_numpy());F[np.ix_(idxs,idxs)].astype('<f8').tofile(w/'stage6_foxp2_gnomad_fin_narrow_ld.bin')
fin.to_csv(o/'stage6_foxp2_gnomad_narrow_variants.csv',index=False);(o/'stage6_foxp2_gnomad_nfe_narrow_ld_receipt.json').write_text(json.dumps({'source':base,'metadata_sha256':hashlib.sha256(mr.content).hexdigest(),'n':n,'blocks':receipt,'whole_matrix_checksum_verified':False},indent=2)+'\n');print('narrow variants',n,flush=True)
