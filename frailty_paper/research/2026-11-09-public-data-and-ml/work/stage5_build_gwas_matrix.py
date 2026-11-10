import csv,collections,pathlib,json,hashlib,datetime
p=pathlib.Path('outputs');queue=list(csv.DictReader(open(p/'stage5_functional_atlas_with_priority_splice.csv')));g=list(csv.DictReader(open(p/'stage5_grch37_gwas_exact_queue.csv')));f=list(csv.DictReader(open(p/'stage5_finngen_insomnia_exact_queue.csv')));credible=list(csv.DictReader(open(p/'stage5_credible_qtl_standardized_deduplicated.csv')))
traits=['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_R13_insomnia'];source={(x['variant_GRCh38'],x['trait_source']):x for x in g};fin={x['variant_GRCh38']:x for x in f};cg=collections.defaultdict(list)
for x in credible:cg[x['variant_GRCh38']].append(x)
fields=['variant_GRCh38_ref_alt','input_classes','input_labels','candidate_gene_context','credible_set_gene_ids','credible_set_tissues','credible_set_record_count','molecular_pip_max_across_sets','avi_phred','ccre_ids']
for t in traits:fields += [t+'_match_status',t+'_beta_alt',t+'_se',t+'_p']
fields+=['interpretation'];rows=[]
for a in queue:
 k=a['variant_GRCh38_ref_alt'];c=cg[k];pip=[float(x['molecular_pip']) for x in c if x['molecular_pip']];r={'variant_GRCh38_ref_alt':k,'input_classes':a['input_classes'],'input_labels':a['input_labels'],'candidate_gene_context':a['candidate_gene_context'],'credible_set_gene_ids':'|'.join(sorted({x['gene_id'] for x in c if x['gene_id']})),'credible_set_tissues':'|'.join(sorted({x['tissue_or_cell'] for x in c if x['tissue_or_cell']})),'credible_set_record_count':len(c),'molecular_pip_max_across_sets':max(pip) if pip else '', 'avi_phred':a['atlas_avi_phred_indexed'],'ccre_ids':a['encode4_ccre_ids']}
 for t in traits:
  x=fin[k] if t=='FinnGen_R13_insomnia' else source[(k,t)];r[t+'_match_status']=x['status'];r[t+'_beta_alt']=x['beta_alt'] if t=='FinnGen_R13_insomnia' else x['gwas_beta_alt'];r[t+'_se']=x['se'] if t=='FinnGen_R13_insomnia' else x['gwas_se'];r[t+'_p']=x['p_value'] if t=='FinnGen_R13_insomnia' else x['gwas_p']
 r['interpretation']='Molecular PIP is not disease PIP; GWAS traits, scales and cohorts differ; exact association does not establish colocalization, common causal variant, or mediation; absent score is not zero.';rows.append(r)
assert len(rows)==len({r['variant_GRCh38_ref_alt'] for r in rows})==488
with open(p/'stage5_variant_multitrait_gwas_matrix.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
stats={t:sum(bool(r[t+'_p']) for r in rows) for t in traits};audit={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'variants':len(rows),'trait_values':stats,'molecular_qtl_variants':sum(bool(cg[r['variant_GRCh38_ref_alt']]) for r in rows),'sha256':hashlib.sha256((p/'stage5_variant_multitrait_gwas_matrix.csv').read_bytes()).hexdigest(),'notes':'Each variant–trait source occurs once; p-values are unadjusted and not evidence of colocalization.'}
(p/'stage5_variant_multitrait_gwas_matrix_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
