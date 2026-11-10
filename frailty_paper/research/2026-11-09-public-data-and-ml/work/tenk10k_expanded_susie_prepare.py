import pandas as pd,numpy as np,json,ast,math
from pathlib import Path
root=Path('.'); base=root/'work/stage6_tenk10k_susie_expanded';base.mkdir(exist_ok=True)
ref=pd.read_csv('outputs/stage6_tenk10k_expanded_eur_ld_order.csv');mat=np.load('outputs/stage6_tenk10k_expanded_eur_ld.npz');
assert np.array_equal(ref.variant_GRCh38.to_numpy(),mat['variant'])
ns={'math':math};tree=ast.parse((root/'work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),'coloc','exec'),ns)
results=[]
for trait,sd in [('Mak_FI',1.1004195728364333),('Mak_FP',.1232461359076205),('FinnGen_insomnia',.2)]:
 d=pd.read_csv(f'outputs/stage6_tenk10k_peak_pal_freq_{trait}_aligned.csv');r=ref.reset_index().merge(d,on='variant_GRCh38',how='inner',validate='one_to_one');r['af_diff']=(r.alt_frequency-r.qtl_alt_af).abs();r=r[r.af_diff<=.1].sort_values('index'); mask=d.variant_GRCh38.isin(r.variant_GRCh38)
 record={'trait':trait,'input_n':len(d),'kept_n':len(r),'reference_af_excluded':int(len(ref.merge(d,on='variant_GRCh38'))-len(r))}
 for name,b,se,prior in [('qtl','qtl_beta_alt','qtl_se',.15),('gwas','gwas_beta_alt','gwas_se',sd)]:
  weights=[ns['bf'](x,y,prior) for x,y in zip(d[b],d[se])];record[name+'_abf_support']=math.exp(ns['logsum']([w for w,m in zip(weights,mask) if m])-ns['logsum'](weights));record[name+'_top_included']=bool(mask.iloc[np.argmax(weights)])
 results.append(record)
 if trait=='Mak_FI':
  r.to_csv(base/'aligned.csv',index=False);idx=r['index'].to_numpy();np.asarray(mat['r'][np.ix_(idx,idx)],dtype='<f8',order='F').ravel(order='F').tofile(base/'ld_f64.bin')
(root/'outputs/stage6_tenk10k_expanded_ld_coverage.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
