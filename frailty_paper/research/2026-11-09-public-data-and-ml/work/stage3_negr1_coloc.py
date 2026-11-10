import csv,math,pathlib,json,statistics,ast,collections
ns={'math':math};p=pathlib.Path
code=ast.parse(p('work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[x for x in code.body if isinstance(x,ast.FunctionDef)],type_ignores=[]),'functions','exec'),ns)
w=[x for x in csv.DictReader(open('work/stage3_negr1_worklist.csv')) if x['candidate']];sm={x['candidate']:x for x in csv.DictReader(open('work/stage3_negr1_finemapping_summary.csv'))};cs=collections.defaultdict(lambda:collections.defaultdict(dict))
for x in csv.DictReader(open('work/stage3_negr1_empty_credible_sets.csv')):
 v=x['variant'];pip=float(x['pip']);old=cs[x['candidate']][x['cs_id']].get(v);assert old is None or abs(pip-old)<1e-12;cs[x['candidate']][x['cs_id']][v]=pip
rows=[];coverage={};skips=[]
for x in w:
 c=x['candidate'];chrom=x['marker_grch38'].split(':')[0].removeprefix('chr');coverage[c]={}
 for trait in ('Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_insomnia'):
  path=p('work/stage3_negr1_qtl/aligned')/(c+'_'+trait+'.csv')
  if not path.exists():skips.append({'candidate':c,'trait':trait,'reason':'no QTL records or matching output'});continue
  d=list(csv.DictReader(path.open()))
  if len(d)<100:skips.append({'candidate':c,'trait':trait,'reason':'fewer than 100 matched SNPs','n':len(d)});continue
  support={'chr'+chrom+'_'+z['position_GRCh38']+'_'+z['ref']+'_'+z['alt'] for z in d};sets={}
  for csid,v in cs[c].items():
   matched=support&v.keys();sets[csid]={'size':len(v),'matched':len(matched),'reported_pip_sum':sum(v.values()),'matched_reported_pip_sum':sum(v[k] for k in matched)}
  if trait=='Atkins_FI':sd=.15*statistics.median(float(z['gwas_se'])*math.sqrt(2*float(z['gwas_N'])*float(z['gwas_af_alt'])*(1-float(z['gwas_af_alt']))) for z in d)
  elif trait=='Jansen_insomnia' or trait=='FinnGen_insomnia':sd=.2
  elif trait=='Mak_FI':sd=1.1004195728364333
  elif trait=='Mak_FP':sd=.1232461359076205
  else:sd=.15*statistics.median(float(z['gwas_se'])*math.sqrt(2*386565*float(z['qtl_maf'])*(1-float(z['qtl_maf']))) for z in d)
  coverage[c][trait]={'n':len(d),'credible_sets':sets,'gwas_prior_sd':sd,'max_abs_gwas_z':max(abs(float(z['gwas_beta_alt'])/float(z['gwas_se'])) for z in d)}
  limitation='no molecular credible set reported' if not sets else 'multiple molecular credible sets; single-signal model inappropriate for decisive interpretation' if len(sets)>1 else 'single molecular set reported; assess matched-variant coverage and GWAS signal multiplicity'
  for scale in (.5,1.,2.):
   l1=[ns['bf'](float(z['qtl_beta_alt']),float(z['qtl_se']),.15) for z in d];l2=[ns['bf'](float(z['gwas_beta_alt']),float(z['gwas_se']),sd*scale) for z in d]
   for p12 in (1e-6,1e-5,5e-5):
    h=ns['posterior'](l1,l2,1e-4,1e-4,p12);rows.append(dict(candidate=c,gene_id=x['gene_id'],gene_symbol=x['gene_symbol'],tissue=x['tissue_v10'],dataset=x['dataset_r7'],trait=trait,n_variants=len(d),qtl_prior_sd=.15,gwas_prior_sd=sd*scale,prior_scale=scale,p1=1e-4,p2=1e-4,p12=p12,H0=h[0],H1=h[1],H2=h[2],H3=h[3],H4=h[4],limitation=limitation))
with open('outputs/stage3_negr1_coloc_sensitivity.csv','w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
p('outputs/stage3_negr1_coloc_coverage.json').write_text(json.dumps({'results':coverage,'skipped':skips},indent=2)+'\n')
print('baseline comparisons',len(rows)//9,'skips',len(skips))
for x in rows:
 if x['prior_scale']==1 and x['p12']==1e-5 and x['trait'] in ('Atkins_FI','Jansen_insomnia'):print(x['candidate'],x['gene_symbol'],x['trait'],round(x['H1'],3),round(x['H3'],3),round(x['H4'],3),x['limitation'])
