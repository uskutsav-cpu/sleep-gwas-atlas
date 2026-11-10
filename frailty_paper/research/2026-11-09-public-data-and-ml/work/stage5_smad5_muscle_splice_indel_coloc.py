import csv,json,math,pathlib,hashlib,ast,collections,gzip
p=pathlib.Path;ns={'math':math};tree=ast.parse(p('work/stage3_smim32_coloc_abf.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),'coloc','exec'),ns)
q=collections.defaultdict(list)
for x in csv.DictReader(open('work/stage3_QTD000285_SMAD5_splice_1mb.tsv'),delimiter='\t'):q[x['variant']].append(x)
c=[x for x in csv.DictReader(open('outputs/stage5_credible_qtl_standardized_deduplicated.csv')) if x['dataset_id']=='QTD000285' and x['gene_id']=='ENSG00000113658'];cs=collections.defaultdict(dict)
for x in c:cs[x['credible_set_id']][x['variant_GRCh38']]=float(x['molecular_pip'])
ind={x['variant_GRCh38']:x for x in c if len(x['effect_allele'])!=len(x['other_allele'])};assert len(ind)==8
match={(x['variant_GRCh38'],x['trait_source']):x for x in csv.DictReader(open('outputs/stage5_indel_grch37_gwas_matches.csv'))};fm={x['variant_GRCh38']:x for x in csv.DictReader(open('outputs/stage5_finngen_insomnia_exact_queue.csv'))};maps={x['variant_GRCh38']:x for x in csv.DictReader(open('outputs/stage5_indel_hg19_mapping.csv'))};raw={}
for t in ['Mak_FI','Mak_FP']:
 wanted={(maps[k]['hg19_chrom'],maps[k]['hg19_position'],frozenset((maps[k]['hg19_ref'],maps[k]['hg19_alt']))) for k in ind}
 with gzip.open('work/stage3_region_cache/'+t+'.regional.tsv.gz','rt') as f:
  for z in csv.DictReader(f,delimiter='\t'):
   key=(z['CHR'],z['POS'],frozenset((z['A1'],z['A2'])))
   if key in wanted:raw[(t,key)]=z
old=list(csv.DictReader(open('outputs/stage3_smad5_splice_coloc_abf_sensitivity.csv')));res=[];coverage=[];exclusions=[]
for trait in ['Mak_FI','Mak_FP','FinnGen_insomnia']:
 src=p('work/stage3_smad5_'+trait+'_aligned.csv');data=list(csv.DictReader(src.open()));baseline=next(x for x in old if x['trait']==trait and float(x['qtl_prior_effect_sd'])==.15 and float(x['p12'])==1e-5);sd=float(baseline['gwas_prior_sd']);prior=ns['posterior']([ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),.15) for x in data],[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x in data],1e-4,1e-4,1e-5);assert max(abs(prior[i]-float(baseline['H'+str(i)])) for i in range(5))<1e-10
 added=[]
 for k,a in ind.items():
  key=k.replace(':','_',2).replace('>','_');v=q[key];assert len({(z['beta'],z['se'],z['pvalue'],z['maf']) for z in v})==1;z=v[0];assert .01<=float(z['maf'])<=.5 and float(z['se'])>0
  if trait=='FinnGen_insomnia':h=fm[k];status=h['status'];beta=h['beta_alt'];se=h['se'];pv=h['p_value'];af=h['af_alt'];ok=status=='exact REF/ALT match'
  else:
   h=match[(k,trait)];status=h['match_status'];beta=h['beta_alt'];se=h['se'];pv=h['p'];ok=status=='unique_exact_allele_match';m=maps[k];mk=(m['hg19_chrom'],m['hg19_position'],frozenset((m['hg19_ref'],m['hg19_alt'])));o=raw.get((trait,mk));af=(float(o['AF1']) if o['A1']==m['hg19_alt'] else 1-float(o['AF1'])) if o else ''
  if not ok or af=='' or not .01<=float(af)<=.99:exclusions.append({'trait':trait,'variant':k,'status':status,'af_alt':af,'reason':'unmatched or outside prior 1% frequency filter'});continue
  item={'variant':z['rsid'],'position_GRCh38':z['position'],'ref':z['ref'],'alt':z['alt'],'qtl_beta_alt':z['beta'],'qtl_se':z['se'],'qtl_p':z['pvalue'],'qtl_maf':z['maf'],'gwas_beta_alt':beta,'gwas_se':se,'gwas_p':pv,'gwas_af_alt':af};added.append(item)
 combined=data+added;keys=[(x['position_GRCh38'],x['ref'],x['alt']) for x in combined];assert len(keys)==len(set(keys));fn=p('outputs/stage5_smad5_splice_'+trait+'_with_credible_indels_aligned.csv')
 with fn.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(combined)
 before={'chr5:'+x['position_GRCh38']+':'+x['ref']+'>'+x['alt'] for x in data};after={'chr5:'+x['position_GRCh38']+':'+x['ref']+'>'+x['alt'] for x in combined};row={'trait':trait,'before_n':len(data),'indels_added':len(added),'after_n':len(combined),'aligned_sha256':hashlib.sha256(fn.read_bytes()).hexdigest()}
 for csid,vals in cs.items():
  label=csid.rsplit('_',1)[-1];row.update({label+'_total':len(vals),label+'_count_before':len(before&vals.keys()),label+'_count_after':len(after&vals.keys()),label+'_pip_before':sum(vals[k] for k in before&vals.keys()),label+'_pip_after':sum(vals[k] for k in after&vals.keys()),label+'_pip_total':sum(vals.values())})
 coverage.append(row)
 for qsd in (.075,.15,.3):
  l1=[ns['bf'](float(x['qtl_beta_alt']),float(x['qtl_se']),qsd) for x in combined];l2=[ns['bf'](float(x['gwas_beta_alt']),float(x['gwas_se']),sd) for x in combined]
  for p12 in (1e-6,1e-5,5e-5):
   h=ns['posterior'](l1,l2,1e-4,1e-4,p12);res.append({'trait':trait,'n':len(combined),'qtl_prior_sd':qsd,'gwas_prior_sd':sd,'p12':p12,**{'H'+str(i):h[i] for i in range(5)}})
for fn,rows in [('stage5_smad5_splice_indel_coloc_sensitivity.csv',res),('stage5_smad5_splice_indel_coloc_coverage.csv',coverage),('stage5_smad5_splice_indel_exclusions.csv',exclusions)]:
 with open('outputs/'+fn,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('coverage',coverage);[print('baseline',x['trait'],x['H1'],x['H3'],x['H4']) for x in res if x['qtl_prior_sd']==.15 and x['p12']==1e-5]
