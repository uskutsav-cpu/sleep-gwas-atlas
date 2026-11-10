import json, math, hashlib, pathlib
import numpy as np
out=[]
for base in ('mak_r13','mak_mvp','mak_ukb'):
 names=('old_FI','old_FFS','mak_FI','mak_FP'); paths=[pathlib.Path('outputs')/f'{base}_{n}.exact.json' for n in names]; d=[json.loads(p.read_text())[0] for p in paths]
 assert len({x['row_order_SNP_sha256'] for x in d})==1
 assert all(x['hsq1']['separators']==x['hsq2']['separators']==d[0]['hsq1']['separators'] for x in d)
 for i,j,label in ((0,1,'original'),(2,3,'Mak')):
  a=np.array(d[i]['delete_rg'],float)[:,0]; b=np.array(d[j]['delete_rg'],float)[:,0]; assert len(a)==len(b)==200
  fac=(len(a)-1)/len(a); va=fac*np.sum((a-a.mean())**2); vb=fac*np.sum((b-b.mean())**2); cv=fac*np.sum((a-a.mean())*(b-b.mean()))
  assert math.isclose(math.sqrt(va),d[i]['rg_se'],rel_tol=1e-10); assert math.isclose(math.sqrt(vb),d[j]['rg_se'],rel_tol=1e-10)
  paired=math.sqrt(fac*np.sum(((b-a)-(b-a).mean())**2)); assert math.isclose(paired**2,va+vb-2*cv,rel_tol=1e-10)
  out.append({'insomnia_source':base,'frailty_pair':label,'genetic_corr_FI':d[i]['rg'],'genetic_corr_physical':d[j]['rg'],'difference_physical_minus_FI':d[j]['rg']-d[i]['rg'],'paired_se':paired,'independence_assumption_se_illustration':math.sqrt(va+vb),'jackknife_corr':cv/math.sqrt(va*vb),'n_snps':d[i]['n_snps'],'input_sha256':[hashlib.sha256(paths[k].read_bytes()).hexdigest() for k in (i,j)]})
p=pathlib.Path('outputs/insomnia_frailty_paired_uncertainty_audit.json');p.write_text(json.dumps({'n_blocks':200,'results':out,'caveat':'Genomic-block jackknife covariance under LDSC assumptions; not a direct measure of participant overlap; naive SE only an illustration'},indent=2)+'\n')
print(json.dumps(out,indent=2))
