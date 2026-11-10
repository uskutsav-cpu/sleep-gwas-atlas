import pathlib,json,csv,math,hashlib
base=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/g/work/sleep-gwas-atlas/frailty_paper/nature_communications_continuation_v6/repro_inputs/native_exact')
paths=[('FinnGen insomnia matched',base/'aux_R13_common_FI_FFS.exact.json'),('MVP insomnia matched',base/'aux_MVP_common_FI_FFS.exact.json'),('FI FFS HFRS common',pathlib.Path('outputs/FI_FFS_HFRS_common.exact.json'))];rows=[]
for label,p in paths:
 d=json.loads(p.read_bytes())[0];a=d['hsq1'];b=d['hsq2'];g=d['gencov'];r=g['total']/math.sqrt(a['total']*b['total']);assert math.isclose(r,d['rg'],rel_tol=1e-12)
 row={'support':label,'n_snps':d['n_snps'],'FI_h2_observed':a['total'],'FI_h2_se':a['total_se'],'FI_h2_precision_z':a['total']/a['total_se'],'FI_h2_intercept':a['intercept'],'FI_intercept_se':a['intercept_se'],'FFS_h2_observed':b['total'],'FFS_h2_se':b['total_se'],'FFS_h2_precision_z':b['total']/b['total_se'],'FFS_h2_intercept':b['intercept'],'FFS_intercept_se':b['intercept_se'],'genetic_covariance':g['total'],'genetic_covariance_se':g['total_se'],'cross_trait_intercept':g['intercept'],'cross_trait_intercept_se':g['intercept_se'],'rg':r,'rg_se':d['rg_se'],'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()};rows.append(row)
with open('outputs/fi_ffs_ldsc_heritability_diagnostics.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
