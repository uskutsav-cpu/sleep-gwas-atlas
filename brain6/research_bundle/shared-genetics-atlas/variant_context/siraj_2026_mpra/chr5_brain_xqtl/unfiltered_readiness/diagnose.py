import csv,json,hashlib,math,numpy as np
from scipy.stats import t,norm
from scipy.optimize import minimize_scalar
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'selected_cpg.tsv').open() as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows)==3335 and all(len(x)==20 and x[10]=='chr5:104484889-104484891:cg24905254' for x in rows)
subset=[]
for x in rows:
 try:
  z,pv=float(x[11]),float(x[4]);beta,se=float(x[12]),float(x[13])
  if se<=0 or not math.isfinite(z) or not (0<pv<=1):continue
  assert math.isclose(beta/se,z,rel_tol=2e-6,abs_tol=2e-7)
  if 2<abs(z)<9:subset.append((z,pv))
 except ValueError:continue
v=np.array(subset);z,pv=v[:,0],v[:,1];assert len(v)==1020
loss=lambda df:float(np.median(abs(np.log(pv)-(np.log(2)+t.logsf(abs(z),df)))))
fit=minimize_scalar(loss,bounds=(50,2000),method='bounded');err=np.log(pv)-(np.log(2)+t.logsf(abs(z),370));assert abs(fit.x-370)<.01 and max(abs(err))<1e-5
with (p.parent/'filtered_region/target_filtered.tsv').open() as f:filtered=list(csv.DictReader(f,delimiter='\t'))
index={(x[3],x[10]):x for x in rows};matches=0
for q in filtered:
 x=index.get((q['variantId'],q['target']))
 if x:
  assert math.isclose(float(x[4]),float(q['pvalue']),rel_tol=1e-12)
  assert math.isclose(float(x[12]),float(q['beta']),rel_tol=1e-12)
  matches+=1
assert matches==121
x=next(x for x in rows if x[3]=='chr5:104568525:A:G');z0=float(x[11]);p0=float(x[4]);res={'source_rows':len(rows),'tested_subset_rule':'2 < abs(beta/SE) < 9 and source p>0','subset_rows':len(subset),'estimated_t_df_descriptive':fit.x,'median_absolute_log_p_error_at_df370':float(np.median(abs(err))),'max_absolute_log_p_error_at_df370':float(max(abs(err))),'rs12658032_source_p':p0,'rs12658032_normal_p':2*norm.sf(abs(z0)),'rs12658032_t_df370_p':2*t.sf(abs(z0),370),'filtered_portal_rows_checked':len(filtered),'filtered_portal_exact_p_beta_matches':matches,'source_output_sha256':hashlib.sha256((p/'selected_cpg.tsv').read_bytes()).hexdigest(),'warning':'df 370 inferred numerically from selected data; source model, independence and full method need independent documentation'}
(p/'diagnosis.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
