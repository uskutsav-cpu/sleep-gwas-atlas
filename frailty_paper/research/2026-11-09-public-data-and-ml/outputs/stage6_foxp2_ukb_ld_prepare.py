#!/usr/bin/env python3
from pathlib import Path
import json,numpy as np,pandas as pd
root=Path(__file__).resolve().parents[1];w=root/'work';o=root/'outputs';meta=json.loads((o/'stage6_foxp2_ukb_ld_metadata.json').read_text());table=meta['table'];ids=[r['SNP'] for r in table];eligible={r['rs']:r for r in json.loads((w/'stage6_foxp2_ukb_eligible.json').read_text())};n=len(ids);assert n==len(eligible)==553
raw=np.fromfile(w/'stage6_foxp2_ukb_ld_f32.bin',dtype='<f4').reshape(n,n).astype(np.float64);sign=np.array([eligible[x]['sign'] for x in ids]);R=raw*sign[:,None]*sign[None,:];assert np.allclose(R,R.T,atol=1e-7) and np.allclose(np.diag(R),1,atol=1e-7)
g0=np.load(w/'stage6_foxp2_1000g_eur_genotypes.npz');lookup={v:i for i,v in enumerate(g0['rsid'])};g=g0['dosage'][[lookup[x] for x in ids]].astype(np.float64);g-=g.mean(axis=1,keepdims=True);g/=np.linalg.norm(g,axis=1,keepdims=True);eur=g@g.T
upper=np.triu_indices(n,1);a=R[upper];b=eur[upper];nz=(abs(a)>.1)&(abs(b)>.1);corr=float(np.corrcoef(a[nz],b[nz])[0,1]);opp=int(((a[nz]*b[nz])<0).sum());e=np.linalg.eigvalsh(R);print('n',n,'min eigen',e.min(),'nonzero sign compare',sum(nz),'opposing',opp,'correlation',corr)
R.astype('<f8').tofile(w/'stage6_foxp2_ukb_ld_float64.bin');data=pd.read_csv(w/'stage6_foxp2_susie_aligned.csv').set_index('rs').loc[ids];data.reset_index().to_csv(w/'stage6_foxp2_ukb_common_aligned.csv',index=False);receipt={'n_snps':n,'minimum_eigenvalue':float(e.min()),'signed_pairs_compare_both_abs_gt_0.1':int(sum(nz)),'opposite_sign_pairs':opp,'signed_r_correlation':corr,'quantization':'source int8, dequantized by LDMatrix','block_diagonal_source':True,'unverified':'Matrix allele convention inferred from metadata A1 and cross-checked against 1000G; study-specific LD and whole archive checksum not available.'};(o/'stage6_foxp2_ukb_ld_qc_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
# Eigenvalue clipping and diagonal normalization provide a PSD sensitivity matrix.
vals,vecs=np.linalg.eigh(R);corrected=(vecs*np.maximum(vals,0))@vecs.T;s=np.sqrt(np.diag(corrected));corrected=corrected/s[:,None]/s[None,:];delta=corrected-R
corrected.astype('<f8').tofile(w/'stage6_foxp2_ukb_ld_psd.bin')
# Exact same support and allele direction for the nested 1000G references.
genotypes=g0['dosage'][[lookup[x] for x in ids]].astype(np.float64);fin_idx=np.loadtxt(w/'stage6_foxp2_fin_reference_indices.txt',dtype=int)-1
for label,subset in [('eur503',genotypes),('fin99',genotypes[:,fin_idx])]:
 centered=subset-subset.mean(axis=1,keepdims=True);norm=np.linalg.norm(centered,axis=1,keepdims=True);assert (norm>0).all();standard=centered/norm;matrix=standard@standard.T;matrix.astype('<f8').tofile(w/f'stage6_foxp2_{label}_ukbcommon_ld.bin')
receipt.update({'negative_eigenvalues':int(sum(vals<0)),'max_absolute_psd_change':float(np.max(abs(delta))),'mean_absolute_psd_change':float(np.mean(abs(delta))),'relative_frobenius_change':float(np.linalg.norm(delta)/np.linalg.norm(R))});(o/'stage6_foxp2_ukb_ld_qc_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
