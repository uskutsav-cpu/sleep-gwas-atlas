import csv,glob,collections,hashlib,pathlib,json,datetime
p=pathlib.Path('outputs');queue=list(csv.DictReader(open(p/'stage5_score_queue.csv')));keys={x['variant_key'] for x in queue}
chroms={'stage3_adarb1_muscle_splice':'21','stage3_foxp2_pituitary':'7','stage3_inka1':'3','stage3_rs323509_testis':'5','stage3_smad5':'5','stage3_smim32':'5','stage3_smim32_muscle':'5','stage4_cd4_smad5':'5'};traits=['Atkins_FI','FinnGen_insomnia','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS'];observations=[];files=[]
for fn in glob.glob('work/*_aligned.csv'):
 stem=pathlib.Path(fn).name
 if '_maf005_' in stem:continue
 match=[(pre,t) for pre in chroms for t in traits if stem==f'{pre}_{t}_aligned.csv']
 if not match:continue
 pre,t=match[0];files.append({'path':fn,'sha256':hashlib.sha256(pathlib.Path(fn).read_bytes()).hexdigest()})
 for x in csv.DictReader(open(fn)):
  key=f"chr{chroms[pre]}:{x['position_GRCh38']}:{x['ref']}>{x['alt']}"
  if key in keys:
   observations.append(dict(variant_GRCh38=key,gwas_trait_source=t,gwas_variant_id=x['variant'],gwas_beta_alt=x['gwas_beta_alt'],gwas_se=x['gwas_se'],gwas_p=x['gwas_p'],gwas_af_alt=x.get('gwas_af_alt',''),gwas_N=x.get('gwas_N',x.get('gwas_n','')),qtl_beta_alt=x.get('qtl_beta_alt',''),qtl_p=x.get('qtl_p',''),region_analysis=pre,source_aligned_file=fn,source_sha256=files[-1]['sha256'],caveat='Filtered and region-specific previously aligned data; source and phenotype differences, cohort overlap and allele/build assumptions apply; not causal or shared-risk evidence'))
assert observations
with open(p/'stage5_credible_variants_gwas_aligned_evidence.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(observations[0]));w.writeheader();w.writerows(observations)
g=collections.defaultdict(list)
for x in observations:g[(x['variant_GRCh38'],x['gwas_trait_source'])].append(x)
conf=[]
for key,rows in g.items():
 vals={(r['gwas_beta_alt'],r['gwas_se'],r['gwas_p']) for r in rows}
 if len(vals)>1:conf.append({'variant':key[0],'trait':key[1],'values':list(vals)})
audit={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_files':len(files),'rows':len(observations),'distinct_variants':len({x['variant_GRCh38'] for x in observations}),'distinct_variant_trait_pairs':len(g),'repeated_rows_across_regional_files':len(observations)-len(g),'conflicting_variant_trait_estimate_groups':conf,'traits':dict(collections.Counter(x['gwas_trait_source'] for x in observations)),'source_files':files,'note':'Exact GRCh38/ref/alt join of pre-existing filtered regional files; absence is not a null result, and duplicate rows are correlated.'}
(p/'stage5_credible_variants_gwas_evidence_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps({k:v for k,v in audit.items() if k not in ['source_files','conflicting_variant_trait_estimate_groups']},indent=2));print('conflict count',len(conf))
