import csv,math,statistics,pathlib,json,hashlib,runpy
f=runpy.run_path('work/stage3_smim32_coloc_abf.py');bf=f['bf'];post=f['posterior'];out=[];diagnostics={}
for src in ('Mak_FI','Mak_FP'):
 p=pathlib.Path('work/stage3_smim32_'+src+'_with_n.csv');rows=list(csv.DictReader(p.open()));assert len(rows)==2404
 def est_sdY(trait):
  xs=[1/float(x[trait+'_se'])**2 for x in rows]
  ys=[]
  for x in rows:
   fr=float(x['qtl_maf']) if trait=='qtl' else float(x['gwas_af_alt']); n=702 if trait=='qtl' else int(x['gwas_n']);ys.append(2*n*fr*(1-fr))
  return math.sqrt(sum(x*y for x,y in zip(xs,ys))/sum(x*x for x in xs))
 sd=est_sdY('gwas');qsd=est_sdY('qtl');diagnostics[src]={'n':len(rows),'estimated_gwas_sdY':sd,'estimated_qtl_sdY':qsd,'n_min':min(int(x['gwas_n']) for x in rows),'n_max':max(int(x['gwas_n']) for x in rows),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 for factor in (.5,1,2):
  for p12 in (1e-6,1e-5,5e-5):
   a=[bf(float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in rows];b=[bf(float(x['gwas_beta_alt']),float(x['gwas_se']),.15*sd*factor) for x in rows];h=post(a,b,1e-4,1e-4,p12)
   out.append({'trait':src,'n_variants':len(rows),'gwas_sdy_estimate':sd,'qtl_effect_prior_sd':.15,'gwas_effect_prior_sd':.15*sd*factor,'gwas_prior_scale_factor':factor,'p1':1e-4,'p2':1e-4,'p12':p12,'H0':h[0],'H1':h[1],'H2':h[2],'H3':h[3],'H4':h[4],'H3_plus_H4':h[3]+h[4],'H4_cond_H3_H4':h[4]/(h[3]+h[4])})
p=pathlib.Path('outputs/stage3_smim32_mak_coloc_abf_sensitivity.csv')
with p.open('w',newline='') as z:w=csv.DictWriter(z,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
pathlib.Path('outputs/stage3_smim32_mak_coloc_receipt.json').write_text(json.dumps({'inputs':diagnostics,'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sdY_estimator':'square root of slope of no-intercept regression 2*N*af*(1-af) on 1/varbeta','limitations':'Estimate residual scale from common variants and sample sizes; beta phenotype transformation not author-confirmed; inferred sample sizes and allele mappings; custom single-signal ABF not validated against coloc; QTL-GWAS overlap unknown; 1Mb restricted coverage; not causal or independent replication'},indent=2)+'\n')
for r in out:
 if r['gwas_prior_scale_factor']==1 and r['p12']==1e-5:print(r)
