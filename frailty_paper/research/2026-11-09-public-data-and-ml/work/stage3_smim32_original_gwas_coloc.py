import csv,math,pathlib,json,statistics,ast
ns={'math':math};p=pathlib.Path
code=ast.parse(p('work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[x for x in code.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),'functions','exec'),ns)
analyses=[('smim32_muscle','outputs/stage3_smim32_muscle_credible_set.csv','5')]
allrows=[];cover={}
for analysis,cs_path,chrom in analyses:
 cs={x['variant'] for x in csv.DictReader(open(cs_path))};cover[analysis]={}
 for trait in ('Atkins_FI','Jansen_insomnia'):
  path=p('work/stage3_'+analysis+'_'+trait+'_aligned.csv');d=list(csv.DictReader(path.open()));support={'chr'+chrom+'_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'] for x in d}
  counts={'n_matched':len(d),'molecular_credible_set_size':len(cs),'molecular_credible_variants_matched':len(cs&support),'max_abs_gwas_z':max(abs(float(x['gwas_beta_alt'])/float(x['gwas_se'])) for x in d)}
  if trait=='Atkins_FI':
   sdY=statistics.median(float(x['gwas_se'])*math.sqrt(2*float(x['gwas_N'])*float(x['gwas_af_alt'])*(1-float(x['gwas_af_alt']))) for x in d);sd=.15*sdY;counts['estimated_sd_from_effective_source_n_and_gwas_frequency']=sdY
  else:sd=.2
  cover[analysis][trait]=counts
  for scale in (.5,1,2):
   l1=[ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in d];l2=[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd*scale) for x in d]
   for p12 in (1e-6,1e-5,5e-5):
    h=ns['posterior'](l1,l2,1e-4,1e-4,p12);allrows.append(dict(analysis=analysis,trait=trait,n_variants=len(d),qtl_prior_sd=.15,gwas_prior_sd=sd*scale,prior_scale=scale,p1=1e-4,p2=1e-4,p12=p12,H0=h[0],H1=h[1],H2=h[2],H3=h[3],H4=h[4]))
out=p('outputs/stage3_smim32_original_gwas_coloc_sensitivity.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
p('outputs/stage3_smim32_original_gwas_coloc_coverage.json').write_text(json.dumps(cover,indent=2)+'\n')
for r in allrows:
 if r['prior_scale']==1 and r['p12']==1e-5:print(r['analysis'],r['trait'],'n',r['n_variants'],'H1',round(r['H1'],4),'H3',round(r['H3'],4),'H4',round(r['H4'],4))
