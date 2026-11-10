import csv,collections,json,pathlib,hashlib,datetime
p=pathlib.Path('outputs');variants=list(csv.DictReader(open(p/'stage5_functional_atlas_with_priority_splice.csv')));vm={x['variant_GRCh38_ref_alt']:x for x in variants};credible=list(csv.DictReader(open(p/'stage5_credible_qtl_standardized_deduplicated.csv')));target=list(csv.DictReader(open(p/'stage5_measured_molecular_evidence.csv')))
by=collections.defaultdict(list)
for r in credible:by[r['variant_GRCh38']].append(('credible',r))
for r in target:
 if r['variant_key'] in vm:by[r['variant_key']].append(('targeted',r))
fields=['variant_GRCh38_ref_alt','input_classes','input_labels','qtl_source_collection','qtl_source','qtl_dataset','gene_id','gene_symbol','tissue_or_cell','molecular_measure','molecular_event','credible_set_id','molecular_pip','effect_allele','other_allele','qtl_alt_effect','qtl_se','qtl_p_value','possible_exact_duplicate_across_collections','direct_variant_disease_association_context','disease_association_scope','avi_phred_indexed','ccre_overlap_ids','ccre_overlap_classes','functional_followup','full_precision_alphagenome_rna_splicing','full_precision_alphagenome_chromatin','prediction_qtl_concordance','interpretation']
rows=[]
for v in variants:
 key=v['variant_GRCh38_ref_alt'];records=by[key] or [('none',{})]
 for kind,x in records:
  cr=kind=='credible';ta=kind=='targeted';event=x.get('molecular_trait_id','') if cr else x.get('molecular_event','');dataset=x.get('dataset_id','') if cr else x.get('dataset','');gene=x.get('gene_id','')
  peers=[y for k,y in records if k!=kind and k!='none' and (y.get('dataset_id') or y.get('dataset'))==dataset and y.get('gene_id')==gene and (y.get('molecular_trait_id') or y.get('molecular_event'))==event] if kind!='none' else []
  # An exact effect is checked before flagging a possible cross-collection duplicate.
  effect=x.get('qtl_beta_alt','') if cr else x.get('effect_alt','');peer_same=False
  for z in peers:
   b=z.get('qtl_beta_alt') or z.get('effect_alt')
   try:peer_same=peer_same or abs(float(effect)-float(b))<1e-8
   except (TypeError,ValueError):pass
  assoc=v['disease_association_context'] if 'exploratory_marker' in v['input_classes'].split('|') else ''
  allele=key.rsplit(':',1)[-1].split('>')
  rows.append(dict(variant_GRCh38_ref_alt=key,input_classes=v['input_classes'],input_labels=v['input_labels'],qtl_source_collection={'credible':'frozen credible-set source','targeted':'targeted QTL inventory','none':'no QTL record in these inventories'}[kind],qtl_source=x.get('source_files','') if cr else x.get('source',''),qtl_dataset=dataset,gene_id=gene,gene_symbol=x.get('gene_symbol',''),tissue_or_cell=x.get('tissue_or_cell',''),molecular_measure=x.get('molecular_measure','') if cr else x.get('molecular_type',''),molecular_event=event,credible_set_id=x.get('credible_set_id',''),molecular_pip=x.get('molecular_pip','') if cr else x.get('credible_set_pip',''),effect_allele=x.get('effect_allele','') if cr else allele[1] if ta else '',other_allele=x.get('other_allele','') if cr else allele[0] if ta else '',qtl_alt_effect=effect,qtl_se=x.get('qtl_se','') if cr else x.get('effect_se',''),qtl_p_value=x.get('qtl_p_value','') if cr else x.get('p_value',''),possible_exact_duplicate_across_collections='yes' if peer_same else 'not identified',direct_variant_disease_association_context=assoc,disease_association_scope='direct marker association as recorded; no QTL colocalization implied' if assoc else 'no direct marker association recorded in this table',avi_phred_indexed=v['atlas_avi_phred_indexed'],ccre_overlap_ids=v['encode4_ccre_ids'],ccre_overlap_classes=v['encode4_ccre_classes'],functional_followup=v['functional_experiment_followup'],full_precision_alphagenome_rna_splicing='not obtained',full_precision_alphagenome_chromatin='not obtained',prediction_qtl_concordance='not assessable: matched full-precision model result unavailable',interpretation='Molecular association and regional context do not establish disease causality; source datasets may overlap; allele-specific experimental validation not established'))
assert set(x['variant_GRCh38_ref_alt'] for x in rows)==set(vm)
with open(p/'stage5_variant_gene_tissue_evidence_atlas.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
counts=collections.Counter(x['qtl_source_collection'] for x in rows);audit={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rows':len(rows),'unique_variants':len(vm),'rows_by_source':counts,'cross_collection_possible_duplicate_rows':sum(x['possible_exact_duplicate_across_collections']=='yes' for x in rows),'variants_with_no_qtl_record_in_either_collection':sum(not by[k] for k in vm),'atlas_sha256':hashlib.sha256((p/'stage5_variant_gene_tissue_evidence_atlas.csv').read_bytes()).hexdigest(),'caveat':'Source QTL rows can be selected, correlated or duplicated across releases; no evidence of model–QTL agreement without matched predictions.'}
(p/'stage5_variant_gene_tissue_evidence_atlas_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
