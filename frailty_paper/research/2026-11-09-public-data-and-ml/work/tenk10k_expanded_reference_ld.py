import csv,collections,json,hashlib
import numpy as np
from pathlib import Path
root=Path(__file__).resolve().parents[1];rows={}
for trait in ['Mak_FI','Mak_FP','FinnGen_insomnia']:
 for r in csv.DictReader((root/f'outputs/stage6_tenk10k_peak_pal_freq_{trait}_aligned.csv').open()):
  c,p,ref,alt=r['variant_GRCh38'].split(':');key=(int(r['gwas_pos37']),frozenset([ref,alt])); rows.setdefault(key,[]).append((r['variant_GRCh38'],ref,alt))
found=collections.defaultdict(list)
with (root/'work/stage4_magma_reference/g1000_eur.bim').open() as f:
 for idx,line in enumerate(f):
  z=line.split()
  if z[0]!='5':continue
  key=(int(z[3]),frozenset((z[4],z[5])))
  if key in rows:found[key].append((idx,z[4],z[5],z[1]))
variants=[]; seen=set()
for key,entries in rows.items():
 unique={x[0] for x in entries}
 if len(unique)!=1 or len(found[key])!=1:continue
 var,ref,alt=entries[0];idx,a1,a2,rsid=found[key][0]
 if var in seen:continue
 seen.add(var);variants.append((var,idx,a1,a2,alt,rsid))
print('union',len(rows),'matched',len(variants),'missing_or_nonunique',len(rows)-len(variants),flush=True)
fam=[x.split()[1] for x in (root/'work/stage4_magma_reference/g1000_eur.fam').open()];N=len(fam);assert N==503 and len(set(fam))==N
step=(N+3)//4;bed=np.memmap(root/'work/stage4_magma_reference/g1000_eur.bed',mode='r',dtype='uint8');assert bytes(bed[:3])==bytes([108,27,1])
indices=np.array([x[1] for x in variants],dtype=np.int64);b=bed[3+indices[:,None]*step+np.arange(step)[None,:]];codes=((b[:,:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(len(variants),-1)[:,:N];G=np.array([2.,np.nan,1.,0.])[codes].T
flip=np.array([x[2]!=x[4] for x in variants]);G[:,flip]=2-G[:,flip];missing=np.isnan(G).sum(axis=0);mean=np.nanmean(G,axis=0);G=np.where(np.isnan(G),mean[None,:],G);sd=np.std(G,axis=0);keep=np.isfinite(sd)&(sd>0);Z=(G[:,keep]-mean[keep])/sd[keep];R=(Z.T@Z/N).astype('float32');v=[x for x,k in zip(variants,keep) if k];assert np.max(np.abs(R-R.T))<1e-6 and np.max(np.abs(np.diag(R)-1))<1e-6
out=root/'outputs/stage6_tenk10k_expanded_eur_ld.npz';np.savez_compressed(out,r=R,variant=np.array([x[0] for x in v]),n_samples=np.int32(N))
with (root/'outputs/stage6_tenk10k_expanded_eur_ld_order.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['variant_GRCh38','bim_index','reference_id','alt_frequency','missing_genotypes']);
 for x,m,miss,k in zip(variants,mean,missing,keep):
  if k:w.writerow([x[0],x[1],x[5],m/2,int(miss)])
rec={'union_keys':len(rows),'unique_matches_before_monomorphic':len(variants),'retained':len(v),'removed_monomorphic':int((~keep).sum()),'reference_n':N,'rank_bound':N-1,'matrix_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'method':'unique b37 position, unordered allele pair, positive-strand liftover from aligned dataset; BIM allele 1 dosage flipped to source ALT; mean imputation; 503 EUR reference','warning':'Palindromic orientation conditional on exact allele and frequency checks; external reference finite sampling; per-marker imputation/source QC incomplete.'}
(root/'outputs/stage6_tenk10k_expanded_eur_ld_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
