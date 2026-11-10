import csv,collections,pathlib,os,json
out=pathlib.Path('outputs');prov=list(csv.DictReader(open(out/'stage5_molecular_credible_set_provenance.csv')));manifest=list(csv.DictReader(open(out/'stage5_variant_input_manifest.csv')));norm={r['source_variant_id']:r for r in csv.DictReader(open(out/'stage5_indel_normalization.csv'))}
source_to_key={}
for x in manifest:
 source=x['source_variant_id']
 if not source:continue
 n=norm.get(source);chrom=x['chr'];pos=n['normalized_position_GRCh38'] if n else x['position_GRCh38'];ref=n['normalized_ref'] if n else x['reference_allele'];alt=n['normalized_alt'] if n else x['alternate_allele'];source_to_key[source]=f'{chrom}:{pos}:{ref}>{alt}'
files=sorted({x['source_file'] for x in prov});rows=[]
for file in files:
 for s in csv.DictReader(open(file)):
  vid=s['variant'];parts=vid.split('_');assert len(parts)==4 and vid in source_to_key
  rows.append({'variant_GRCh38_normalized':source_to_key[vid],'source_variant_id':vid,'source_reference':parts[2],'source_alternate':parts[3],'gene_id':s.get('gene_id',''),'gene_symbol_source':s.get('gene_symbol',''),'molecular_trait_id':s.get('molecular_trait_id',''),'credible_set_id':s.get('cs_id',''),'dataset_id':s.get('dataset',''),'source_phenotype_label':s.get('qtl_phenotype',''),'candidate':s.get('candidate',''),'molecular_pip':s.get('pip',''),'qtl_p_value':s.get('pvalue',''),'qtl_beta_source':s.get('beta',''),'qtl_se_source':s.get('se',''),'qtl_z_source':s.get('z',''),'source_file':file,'effect_orientation_status':'source alleles preserved; confirm dataset-specific effect allele before directional model comparison','caveat':'fine-mapped molecular QTL association, not phenotype causality; same variant may occur across multiple sets and related datasets'})
assert len(rows)==len(prov)==651
with open(out/'stage5_credible_set_qtl_effects.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
p=out/'stage5_functional_atlas_with_priority_splice.csv';a=list(csv.DictReader(open(p)));counts=collections.Counter(r['variant_GRCh38_normalized'] for r in rows);assert set(counts).issubset({x['variant_GRCh38_ref_alt'] for x in a});fields=list(a[0])+['frozen_credible_set_qtl_effect_record_count']
for x in a:x['frozen_credible_set_qtl_effect_record_count']=counts[x['variant_GRCh38_ref_alt']]
with open(str(p)+'.tmp','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(a)
os.replace(str(p)+'.tmp',p)
print('Saved',len(rows),'source records across',len(counts),'variants; existing measured evidence table preserved')
