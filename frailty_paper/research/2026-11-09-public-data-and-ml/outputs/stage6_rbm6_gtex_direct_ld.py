from pathlib import Path
import sys,json,hashlib,requests,numpy as np,pandas as pd
O=Path('outputs');W=Path('work');sys.path.insert(0,str(W/'ldcov'))
from ldcov.io.blockmatrix.partitioner import GridPartitioner
from ldcov.io.blockmatrix.block_codec import decode_block
pop=sys.argv[1];assert pop in ('fin','nfe');d=pd.read_csv(O/'stage6_rbm6_gtex_direct_aligned.csv');base=f'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.adj.ld.bm';meta_r=requests.get(base+'/metadata.json',timeout=45);meta_r.raise_for_status();m=meta_r.json();bs=m['blockSize'];p=GridPartitioner(bs,m['nRows'],m['nCols'],m.get('maybeFiltered'));ix=d[pop+'_idx'].to_numpy();sgn=d[pop+'_sign'].to_numpy();bid=ix//bs;groups=sorted(set(bid));M=np.full((len(d),len(d)),np.nan);receipt=[]
for ai,a in enumerate(groups):
 for b in groups[ai:]:
  slot=p.part_slot(int(a),int(b));
  if slot is None:raise ValueError('reference matrix block absent; do not infer zeros')
  name=m['partFiles'][slot];r=requests.get(base+'/parts/'+name,timeout=180);r.raise_for_status();raw=r.content;B=decode_block(raw);i=np.flatnonzero(bid==a);j=np.flatnonzero(bid==b);sub=B[np.ix_(ix[i]%bs,ix[j]%bs)].copy()
  if a==b:sub=sub+sub.T;np.fill_diagonal(sub,np.diag(B)[ix[i]%bs])
  M[np.ix_(i,j)]=sub
  if a!=b:M[np.ix_(j,i)]=sub.T
  receipt.append({'block':[int(a),int(b)],'name':name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()});print(pop,'read',a,b,flush=True);del r,raw,B,sub
M*=sgn[:,None]*sgn[None,:];assert np.isfinite(M).all() and np.allclose(M,M.T,atol=1e-7) and np.allclose(np.diag(M),1,atol=1e-7);eig=np.linalg.eigvalsh(M);print('min eig',eig[0],flush=True)
if eig[0]<-1e-8:raise ValueError('reference matrix not PSD; do not apply automatic correction')
M.astype('<f8').tofile(W/f'stage6_rbm6_gtex_direct_{pop}_ld.bin');(O/f'stage6_rbm6_gtex_direct_{pop}_ld_receipt.json').write_text(json.dumps({'source':base,'metadata_sha256':hashlib.sha256(meta_r.content).hexdigest(),'n':len(d),'blocks':receipt,'min_eigenvalue':float(eig[0]),'full_remote_matrix_checksum_verified':False},indent=2)+'\n')
