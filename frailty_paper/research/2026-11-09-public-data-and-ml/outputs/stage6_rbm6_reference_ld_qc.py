#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,json
root=Path(__file__).resolve().parents[1];out=root/'outputs';work=root/'work';d=np.load(out/'stage6_rbm6_1000g_eur_reference_ld.npz');ids=d['rsid'];r=d['ld'].astype(np.float64);rows=pd.read_csv(out/'stage6_rbm6_five_gwas_common_region.csv.gz').set_index('rsid').loc[ids];assert r.shape==(len(rows),len(rows));i,j=np.where(np.triu(np.abs(r)>.98,1));report=[];pairs=[]
for name in ['Jansen_insomnia','Mak_FI','FinnGen_HFRS']:
 z=(rows[name+'_beta_alt']/rows[name+'_se']).to_numpy();delta=np.abs(z[i]-np.sign(r[i,j])*z[j]);report.append({'trait':name,'n_snps':len(ids),'high_ld_pairs':len(i),'pairs_abs_z_diff_gt_2':int(sum(delta>2)),'pairs_abs_z_diff_gt_3':int(sum(delta>3)),'pairs_abs_z_diff_gt_5':int(sum(delta>5)),'max_abs_z_diff':float(delta.max()) if len(delta) else None});
 for k in np.argsort(delta)[-20:][::-1]:pairs.append({'trait':name,'rsid1':ids[i[k]],'rsid2':ids[j[k]],'reference_r':r[i[k],j[k]],'z1':z[i[k]],'z2':z[j[k]],'aligned_abs_z_difference':delta[k]})
pd.DataFrame(pairs).to_csv(out/'stage6_rbm6_high_ld_z_mismatch.csv',index=False);(out/'stage6_rbm6_high_ld_z_diagnostics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
# Export row-major covariance; symmetry makes it identical in R's column-major read.
r.tofile(work/'stage6_rbm6_eur_ld_float64.bin');rows.reset_index().to_csv(work/'stage6_rbm6_ld_gwas_aligned.csv',index=False)
