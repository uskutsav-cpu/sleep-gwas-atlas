import csv,json,math,pathlib,hashlib,ast,collections,gzip
p=pathlib.Path;ns={'math':math};tree=ast.parse(p('work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'coloc','exec'),ns)
q=collections.defaultdict(list)
for x in csv.DictReader(open('work/stage4_QTD000612_SMAD5_1mb.tsv'),delimiter='\t'):q[x['variant']].append(x)
cs={}
for x in csv.DictReader(open('outputs/stage4_cd4_smad5_credible_set.csv')):
 assert x['variant'] not in cs or cs[x['variant']]==float(x['pip']);cs[x['variant']]=float(x['pip'])
ind=[x for x in csv.DictReader(open('outputs/stage5_credible_qtl_standardized_deduplicated.csv')) if x['dataset_id']=='QTD000612' and len(x['effect_allele'])!=len(x['other_allele'])];assert len(ind)==11
match={(x['variant_GRCh38'],x['trait_source']):x for x in csv.DictReader(open('outputs/stage5_indel_grch37_gwas_matches.csv'))};fm={x['variant_GRCh38']:x for x in csv.DictReader(open('outputs/stage5_finngen_insomnia_exact_queue.csv'))};maps={x['variant_GRCh38']:x for x in csv.DictReader(open('outputs/stage5_indel_hg19_mapping.csv'))}
raw={}
for t in ['Mak_FI','Mak_FP']:
 wanted={(maps[x['variant_GRCh38']]['hg19_chrom'],maps[x['variant_GRCh38']]['hg19_position'],frozenset((maps[x['variant_GRCh38']]['hg19_ref'],maps[x['variant_GRCh38']]['hg19_alt']))) for x in ind}
 with gzip.open('work/stage3_region_cache/'+t+'.regional.tsv.gz','rt') as f:
  for z in csv.DictReader(f,delimiter='\t'):
   key=(z['CHR'],z['POS'],frozenset((z['A1'],z['A2'])))
   if key in wanted:raw[(t,key)]=z
old=list(csv.DictReader(open('outputs/stage4_cd4_smad5_coloc_abf_sensitivity.csv')));results=[];coverage=[];exclusions=[]
for trait,sd in [('Mak_FI',1.1004195728364333),('Mak_FP',.1232461359076205),('FinnGen_insomnia',.2)]:
 src=p('work/stage4_cd4_smad5_'+trait+'_aligned.csv');data=list(csv.DictReader(src.open()));good=[]
 for x in data:
  k='chr5_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'];values=q.get(k,[]);r=float(values[0]['r2']) if values and values[0]['r2']!='NA' else float('nan')
  if math.isfinite(r) and r>=.8:good.append(x)
 prior=next(x for x in old if x['trait']==trait and float(x['qtl_prior_sd'])==.15 and float(x['p12'])==1e-5);prev=ns['posterior']([ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in good],[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x in good],1e-4,1e-4,1e-5);assert max(abs(prev[i]-float(prior['H'+str(i)])) for i in range(5))<1e-10
 added=[]
 for a in ind:
  k=a['variant_GRCh38'];source=k.replace(':','_',2).replace('>','_');vals=q[source];assert len({(z['beta'],z['se'],z['pvalue'],z['maf'],z['r2']) for z in vals})==1;z=vals[0];assert float(z['r2'])>=.8 and .01<=float(z['maf'])<=.5
  if trait=='FinnGen_insomnia':
   h=fm[k];status=h['status'];beta=h['beta_alt'];se=h['se'];pv=h['p_value'];af=h['af_alt'];ok=status=='exact REF/ALT match'
  else:
   h=match[(k,trait)];status=h['match_status'];beta=h['beta_alt'];se=h['se'];pv=h['p'];ok=status=='unique_exact_allele_match';m=maps[k];mk=(m['hg19_chrom'],m['hg19_position'],frozenset((m['hg19_ref'],m['hg19_alt'])));o=raw.get((trait,mk));af=(float(o['AF1']) if o['A1']==m['hg19_alt'] else 1-float(o['AF1'])) if o else ''
  if not ok or not af or not .01<=float(af)<=.99:exclusions.append({'trait':trait,'variant':k,'status':status,'af':af});continue
  item={'variant':z['rsid'],'position_GRCh38':z['position'],'ref':z['ref'],'alt':z['alt'],'qtl_beta_alt':z['beta'],'qtl_se':z['se'],'qtl_p':z['pvalue'],'qtl_maf':z['maf'],'gwas_beta_alt':beta,'gwas_se':se,'gwas_p':pv,'gwas_af_alt':af};added.append(item)
 combined=good+added;keys=[(x['position_GRCh38'],x['ref'],x['alt']) for x in combined];assert len(keys)==len(set(keys));fn=p('outputs/stage5_cd4_smad5_'+trait+'_with_credible_indels_aligned.csv')
 with fn.open('w',newline='') as z:w=csv.DictWriter(z,fieldnames=list(good[0]));w.writeheader();w.writerows(combined)
 before={'chr5_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'] for x in good};after={'chr5_'+x['position_GRCh38']+'_'+x['ref']+'_'+x['alt'] for x in combined};coverage.append({'trait':trait,'baseline_n':len(good),'n_indels_added':len(added),'updated_n':len(combined),'credible_count_before':len(before&cs.keys()),'credible_count_after':len(after&cs.keys()),'credible_pip_before':sum(cs[x] for x in before&cs.keys()),'credible_pip_after':sum(cs[x] for x in after&cs.keys()),'total_credible_pip':sum(cs.values()),'aligned_input_sha256':hashlib.sha256(fn.read_bytes()).hexdigest()})
 for sdq in (.075,.15,.30):
  l1=[ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),sdq) for x in combined];l2=[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x in combined]
  for p12 in (1e-6,1e-5,5e-5):
   h=ns['posterior'](l1,l2,1e-4,1e-4,p12);results.append({'trait':trait,'n':len(combined),'qtl_prior_sd':sdq,'gwas_prior_sd':sd,'p12':p12,**{'H'+str(i):h[i] for i in range(5)}})
for file,vals in [('stage5_cd4_smad5_indel_coloc_sensitivity.csv',results),('stage5_cd4_smad5_indel_coloc_coverage.csv',coverage),('stage5_cd4_smad5_indel_exclusions.csv',exclusions)]:
 with open('outputs/'+file,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(vals[0]));w.writeheader();w.writerows(vals)
print('coverage',coverage);print('baseline');[print(x['trait'],x['H1'],x['H3'],x['H4']) for x in results if x['qtl_prior_sd']==.15 and x['p12']==1e-5]
