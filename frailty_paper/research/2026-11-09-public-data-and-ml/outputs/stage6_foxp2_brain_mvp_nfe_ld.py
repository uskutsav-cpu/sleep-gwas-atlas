from pathlib import Path
import sys,json,hashlib,requests,numpy as np,pandas as pd
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';sys.path.insert(0,str(w/'ldcov'))
from ldcov.io.blockmatrix.partitioner import GridPartitioner
from ldcov.io.blockmatrix.block_codec import decode_block
idx=pd.read_csv(o/'stage6_foxp2_brain_mvp_nfe_variant_index.csv');src=pd.read_csv(w/'stage6_foxp2_brain_mvp_pair_input.csv');d=idx.merge(src,on='rs',validate='one_to_one',suffixes=('','_source'));d['reference_anchor_af']=np.where(d.sign_to_anchor==1,d.gnomad_af_alt,1-d.gnomad_af_alt);d['af_difference']=abs(d.reference_anchor_af-d.effect_allele_frequency);d=d[d.af_difference<=.2].sort_values('gnomad_index').reset_index(drop=True);n=len(d)
base='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.bm';r=requests.get(base+'/metadata.json',timeout=40);r.raise_for_status();meta=r.json();bs=meta['blockSize'];p=GridPartitioner(bs,meta['nRows'],meta['nCols'],meta.get('maybeFiltered'));ix=d.gnomad_index.to_numpy();bid=ix//bs;groups=sorted(set(bid));M=np.full((n,n),np.nan);receipt=[]
for aa,a in enumerate(groups):
 for b in groups[aa:]:
  slot=p.part_slot(int(a),int(b));
  if slot is None:raise ValueError('Missing matrix block')
  name=meta['partFiles'][slot];r=requests.get(base+'/parts/'+name,timeout=180);r.raise_for_status();raw=r.content;B=decode_block(raw);i=np.flatnonzero(bid==a);j=np.flatnonzero(bid==b);sub=B[np.ix_(ix[i]%bs,ix[j]%bs)].copy()
  if a==b:sub=sub+sub.T;np.fill_diagonal(sub,np.diag(B)[ix[i]%bs])
  M[np.ix_(i,j)]=sub
  if a!=b:M[np.ix_(j,i)]=sub.T
  receipt.append({'block':[int(a),int(b)],'name':name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()});print('read block',a,b,flush=True);del B,raw,r,sub
sgn=d.sign_to_anchor.to_numpy();M*=sgn[:,None]*sgn[None,:];assert np.isfinite(M).all() and np.allclose(M,M.T,atol=1e-7) and np.allclose(np.diag(M),1,atol=1e-7);ev=np.linalg.eigvalsh(M);assert ev[0]>-1e-8;M.astype('<f8').tofile(w/'stage6_foxp2_brain_mvp_nfe_ld.bin');d.to_csv(o/'stage6_foxp2_brain_mvp_nfe_aligned.csv',index=False);(o/'stage6_foxp2_brain_mvp_nfe_ld_receipt.json').write_text(json.dumps({'source':base,'n':n,'AF_excluded':len(idx)-n,'min_eigenvalue':float(ev[0]),'AF_difference_max':float(d.af_difference.max()),'blocks':receipt,'whole_matrix_checksum_verified':False},indent=2)+'\n');print('variants',n,'min eigen',ev[0],flush=True)
