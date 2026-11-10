import csv, json, math, pathlib, hashlib, ast
D=pathlib.Path
ns={'math':math}
tree=ast.parse(D('work/stage3_smim32_coloc_abf.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'coloc_functions','exec'),ns)
q={}
for x in csv.DictReader(open('work/stage4_QTD000612_SMAD5_1mb.tsv'),delimiter='\t'):
 k=x['variant']; r=float(x['r2']) if x['r2']!='NA' else float('nan')
 if k in q: assert q[k]==r
 q[k]=r
cs={}
for x in csv.DictReader(open('outputs/stage4_cd4_smad5_credible_set.csv')):
 k=x['variant'];p=float(x['pip']);assert k not in cs or cs[k]==p;cs[k]=p
assert len(cs)==70
rows=[]; summary={}
for trait,sd in [('Mak_FI',1.1004195728364333),('Mak_FP',0.1232461359076205),('FinnGen_insomnia',.2)]:
 path=D('work/stage4_cd4_smad5_'+trait+'_aligned.csv');data=list(csv.DictReader(path.open()));good=[]
 for x in data:
  k='chr5_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'];r=q.get(k,float('nan'))
  if math.isfinite(r) and r>=.8:good.append((x,k))
 covered={k for x,k in good}&cs.keys(); raw={ 'chr5_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'] for x in data}&cs.keys()
 summary[trait]={'input_variants':len(data),'qtl_r2_ge_0_8_variants':len(good),'credible_set_variants_covered_before_r2':len(raw),'credible_set_variants_covered_after_r2':len(covered),'credible_set_pip_covered_after_r2':sum(cs[k] for k in covered),'credible_set_total_pip':sum(cs.values()),'max_abs_gwas_z':max(abs(float(x['gwas_beta_alt'])/float(x['gwas_se'])) for x,k in good),'candidate_retained':any(x['variant']=='rs7734045' for x,k in good),'aligned_input_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
 for qsd in (.075,.15,.3):
  l1=[ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),qsd) for x,k in good];l2=[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x,k in good]
  for p12 in (1e-6,1e-5,5e-5):
   h=ns['posterior'](l1,l2,1e-4,1e-4,p12)
   rows.append(dict(trait=trait,n_variants=len(good),qtl_prior_sd=qsd,gwas_prior_sd=sd,p1=1e-4,p2=1e-4,p12=p12,H0=h[0],H1=h[1],H2=h[2],H3=h[3],H4=h[4],H4_cond_H3_H4=h[4]/(h[3]+h[4])))
out=D('outputs/stage4_cd4_smad5_coloc_abf_sensitivity.csv')
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
D('outputs/stage4_cd4_smad5_coloc_coverage.json').write_text(json.dumps({'traits':summary,'qtl_credible_set_variants':len(cs),'method':'single-signal ABF, non-palindromic biallelic exact matched SNPs with MAF >= .01 and QTL imputation r2 >= .8','window':'chr5:135702653-136702653 GRCh38','warning':'exploratory and prior sensitive; association does not establish mediation or shared insomnia-frailty risk; disease-signal multiplicity, sample overlap and coverage not fully established'},indent=2)+'\n')
print(json.dumps(summary,indent=2));print(*[x for x in rows if x['p12']==1e-5 and x['qtl_prior_sd']==.15],sep='\n')
