import csv,collections,json,pathlib,datetime
p=pathlib.Path('outputs');r=list(csv.DictReader(open(p/'stage5_credible_set_qtl_effects.csv')))
tissues={'QTD000126':('Adrenal_Gland','gene expression'),'QTD000141':('Artery_Tibial','gene expression'),'QTD000151':('Brain_Anterior_cingulate_cortex_BA24','gene expression'),'QTD000191':('Brain_Nucleus_accumbens_basal_ganglia','gene expression'),'QTD000281':('Muscle_Skeletal','gene expression'),'QTD000285':('Muscle_Skeletal','LeafCutter splicing'),'QTD000301':('Pituitary','gene expression'),'QTD000340':('Testis','LeafCutter splicing'),'QTD000612':('CD4_Naive','pseudobulk gene expression')}
def dataset(x):
 f=x['source_file'];phen=x['source_phenotype_label']
 if f.endswith('stage3_skeletal_qtl_credible_set_variants.csv'):return 'QTD000281' if phen=='SMIM32 expression' else 'QTD000285'
 if f.endswith('stage3_smim32_muscle_credible_set.csv') or f.endswith('stage3_inka1_muscle_credible_set.csv'):return 'QTD000281'
 if f.endswith('stage3_adarb1_muscle_splice_credible_set.csv'):return 'QTD000285'
 if f.endswith('stage3_foxp2_pituitary_susie_credible_set.csv'):return 'QTD000301'
 if f.endswith('stage3_rs323509_testis_splice_credible_set.csv'):return 'QTD000340'
 if f.endswith('stage4_cd4_smad5_credible_set.csv'):return 'QTD000612'
 return x['dataset_id']
g=collections.defaultdict(list)
for x in r:g[(dataset(x),x['molecular_trait_id'],x['credible_set_id'],x['source_variant_id'])].append(x)
rows=[]
for (d,event,cs,vid),same in g.items():
 a=same[0];assert d in tissues and len({(x['qtl_beta_source'],x['qtl_se_source'],x['qtl_p_value'],x['molecular_pip']) for x in same})==1
 chrom,pos,ref,alt=vid.split('_');assert a['variant_GRCh38_normalized']==f'{chrom}:{pos}:{ref}>{alt}'
 rows.append(dict(variant_GRCh38=a['variant_GRCh38_normalized'],dataset_id=d,tissue_or_cell=tissues[d][0],molecular_measure=tissues[d][1],gene_id=a['gene_id'],molecular_trait_id=event,credible_set_id=cs,molecular_pip=a['molecular_pip'],effect_allele=alt,other_allele=ref,qtl_beta_alt=a['qtl_beta_source'],qtl_se=a['qtl_se_source'],qtl_p_value=a['qtl_p_value'],source_files='|'.join(sorted({x['source_file'] for x in same})),source_row_copies=len(same),allele_orientation_basis='eQTL Catalogue ALT effect allele; source and normalized representations identical',caveat='Molecular QTL, not disease causal probability; overlapping samples and model training possible; splicing effects are event-specific'))
assert len(r)==651 and len(rows)==525 and len({x['variant_GRCh38'] for x in rows})==444
with open(p/'stage5_credible_qtl_standardized_deduplicated.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
a=list(csv.DictReader(open(p/'stage5_measured_molecular_evidence.csv'))); comparisons=[]
for x in rows:
 for y in a:
  if y['variant_key']==x['variant_GRCh38'] and y['gene_id']==x['gene_id'] and y['molecular_event']==x['molecular_trait_id'] and y['dataset']==x['dataset_id'] and y['effect_alt']:
   comparisons.append(abs(float(y['effect_alt'])-float(x['qtl_beta_alt'])))
assert len(comparisons)>=1 and max(comparisons)<1e-8
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_rows':len(r),'deduplicated_rows':len(rows),'duplicate_export_rows':len(r)-len(rows),'distinct_variants':len({x['variant_GRCh38'] for x in rows}),'duplicate_conflicts':0,'exact_crosscheck_count':len(comparisons),'max_beta_absolute_difference':max(comparisons),'source':'https://www.ebi.ac.uk/eqtl/Data_access/','caveat':'Coverage and allele bookkeeping checked; tissue labels inherited from release metadata, biological comparability and independence not established.'}
(p/'stage5_credible_qtl_standardization_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
