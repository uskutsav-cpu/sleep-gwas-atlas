import csv,gzip,hashlib,json,pathlib,collections,datetime,os
p=pathlib.Path('outputs');q=list(csv.DictReader(open(p/'stage5_score_queue.csv')));keys={x['variant_key'] for x in q};cache=pathlib.Path('work/stage3_region_cache/FinnGen_insomnia.regional.tsv.gz');found={}
with gzip.open(cache,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  k=f"chr{x['#chrom']}:{x['pos']}:{x['ref']}>{x['alt']}"
  if k in keys:
   if k in found:raise ValueError('duplicate exact allele '+k)
   found[k]=x
previous=[x for x in csv.DictReader(open(p/'stage5_credible_variants_gwas_aligned_evidence.csv')) if x['gwas_trait_source']=='FinnGen_insomnia'];differences=[]
for x in previous:
 a=found.get(x['variant_GRCh38'])
 if not a: differences.append((x['variant_GRCh38'],'missing'))
 elif abs(float(a['beta'])-float(x['gwas_beta_alt']))>1e-8 or abs(float(a['pval'])-float(x['gwas_p']))>1e-8:differences.append((x['variant_GRCh38'],'different'))
rows=[]
for x in q:
 k=x['variant_key'];z=found.get(k)
 rows.append(dict(variant_GRCh38=k,input_classes=x['input_classes'],source='FinnGen R13 F5_INSOMNIA regional cache',status='exact REF/ALT match' if z else 'not present in selected regional cache',rsids=z['rsids'] if z else '',beta_alt=z['beta'] if z else '',se=z['sebeta'] if z else '',p_value=z['pval'] if z else '',af_alt=z['af_alt'] if z else '',af_alt_cases=z['af_alt_cases'] if z else '',af_alt_controls=z['af_alt_controls'] if z else '',nearest_genes_source=z['nearest_genes'] if z else '',caveat='Disease association, not causal or shared-risk proof; absent means no exact row in cache; original source reference/build presumed GRCh38'))
with open(p/'stage5_finngen_insomnia_exact_queue.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
a=list(csv.DictReader(open(p/'stage5_variant_gene_tissue_evidence_atlas.csv')));fields=list(a[0])+['finngen_r13_insomnia_exact_p','finngen_r13_insomnia_beta_alt','finngen_r13_insomnia_match_status']
for x in a:
 z=found.get(x['variant_GRCh38_ref_alt']);x['finngen_r13_insomnia_exact_p']=z['pval'] if z else '';x['finngen_r13_insomnia_beta_alt']=z['beta'] if z else '';x['finngen_r13_insomnia_match_status']='exact REF/ALT match' if z else 'no exact record in regional cache'
path=p/'stage5_variant_gene_tissue_evidence_atlas.csv'
with open(str(path)+'.tmp','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(a)
os.replace(str(path)+'.tmp',path)
receipt={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'queued_variants':len(q),'exact_matches':len(found),'not_found_in_cache':len(q)-len(found),'indel_matches':sum(k in found and len(x['reference_allele'])!=len(x['alternate_allele']) for x in q for k in [x['variant_key']]),'previous_aligned_rows_checked':len(previous),'discrepancies':differences,'cache_sha256':hashlib.sha256(cache.read_bytes()).hexdigest(),'queue_file':str(p/'stage5_score_queue.csv'),'source_metadata':'stage3_gwas_regional_cache_receipt.json','notes':'Exact REF/ALT matching does not independently verify original source build or biological causality.'}
(p/'stage5_finngen_insomnia_exact_queue_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');audit=p/'stage5_variant_gene_tissue_evidence_atlas_audit.json';d=json.loads(audit.read_text());d['created_utc']=receipt['date_utc'];d['atlas_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();d['finngen_exact_variant_matches']=len(found);audit.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(receipt,indent=2));assert len(q)==488 and len(found)==474 and not differences
