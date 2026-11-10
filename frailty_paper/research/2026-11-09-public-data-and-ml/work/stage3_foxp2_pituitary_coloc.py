import csv,math,pathlib,json,statistics,ast
ns={'math':math};p=pathlib.Path
code=ast.parse(p('work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[x for x in code.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),'functions','exec'),ns)
cs={x['variant'] for x in csv.DictReader(open('outputs/stage3_foxp2_pituitary_susie_credible_set.csv'))};assert len(cs)==4
rows=[];cover={}
for trait in ('Mak_FI','Mak_FP','Ye_FFS','FinnGen_insomnia'):
 d=list(csv.DictReader(open('work/stage3_foxp2_pituitary_'+trait+'_aligned.csv')));support={'chr7_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'] for x in d};cover[trait]={'n_matched':len(d),'molecular_cs_variants':len(cs),'molecular_cs_variants_matched':len(cs&support),'candidate_retained':any(x['variant']in ('rs12536335','rs12671330','rs17312861') for x in d),'max_abs_gwas_z':max(abs(float(x['gwas_beta_alt'])/float(x['gwas_se'])) for x in d)}
 if trait=='Ye_FFS':
  estimates=[float(x['gwas_se'])*math.sqrt(2*386565*float(x['qtl_maf'])*(1-float(x['qtl_maf']))) for x in d]
  sd=.15*statistics.median(estimates);cover[trait]['estimated_phenotype_sd_qtl_maf_proxy']=statistics.median(estimates)
 elif trait=='Mak_FI':sd=1.1004195728364333
 elif trait=='Mak_FP':sd=.1232461359076205
 else:sd=.2
 for scale in (.5,1.,2.):
  l1=[ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in d];l2=[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd*scale) for x in d]
  for p12 in (1e-6,1e-5,5e-5):
   h=ns['posterior'](l1,l2,1e-4,1e-4,p12);rows.append(dict(trait=trait,n_variants=len(d),qtl_prior_sd=.15,gwas_prior_sd=sd*scale,prior_scale=scale,p1=1e-4,p2=1e-4,p12=p12,H0=h[0],H1=h[1],H2=h[2],H3=h[3],H4=h[4]))
with open('outputs/stage3_foxp2_pituitary_coloc_sensitivity.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
p('outputs/stage3_foxp2_pituitary_coloc_coverage.json').write_text(json.dumps(cover,indent=2)+'\n')
print(json.dumps(cover,indent=2));print(*[x for x in rows if x['p12']==1e-5 and x['prior_scale']==1],sep='\n')
