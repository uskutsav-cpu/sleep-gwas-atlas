import csv,collections
edges=list(csv.DictReader(open('outputs/stage4_variant_gene_phenotype_network_scaffold.csv')));fields=list(edges[0])
def add(a,b,typ,source,trait,tissue,p,details):edges.append(dict(zip(fields,[a,b,typ,source,trait,tissue,p,details])))
for x in csv.DictReader(open('outputs/stage5_preliminary_functional_atlas.csv')):
 allele='allele_GRCh38:'+x['variant_GRCh38_ref_alt']
 for label in x['input_labels'].split('|'):
  if label.startswith('rs'):add('variant:'+label,allele,'verified_candidate_allele_mapping','Ensembl GRCh37→38 and GRCh38 reference','allelic mapping','','','reference/alternate pair and liftover checked; QTL or phenotype causality not implied')
 if x['atlas_avi_phred_displayed']!='':add(allele,'prediction:AlphaGenome_Atlas_AVI','AlphaGenome_Atlas_website_AVI','Google DeepMind AlphaGenome Atlas web view','broad variant-impact score','','','displayed rounded Phred '+x['atlas_avi_phred_displayed']+'; top genome-wide SNV percent '+x['atlas_top_genomewide_percent_displayed']+'; '+x['atlas_web_url']+'; not a phenotype risk or validation')
for x in csv.DictReader(open('outputs/stage5_initial_atlas_tissue_observations.csv')):
 allele='allele_GRCh38:'+x['variant_id'];t='cell_or_tissue:AlphaGenome:'+x['biosample'];sc=x['scorer'];detail='displayed raw '+x['displayed_raw_score']+'; website quantile display '+x['displayed_website_quantile']+'; assay '+x['assay']+'; source '+x['source_url']+'; '+x['caveat']
 if x['reported_gene'] in ('FOXP2','SMAD5'):
  gid={'FOXP2':'ENSG00000128573','SMAD5':'ENSG00000113658'}[x['reported_gene']];add(allele,'gene:'+gid,'AlphaGenome_website_gene_expression_score','Google DeepMind AlphaGenome Atlas',sc,x['biosample'],'',detail)
 else:add(allele,t,'AlphaGenome_website_local_regulatory_score','Google DeepMind AlphaGenome Atlas',sc,x['biosample'],'',detail)
for x in csv.DictReader(open('outputs/stage5_rs323509_atlas_initial_observations.csv')):
 if x['output_type']!='RNA-seq':continue
 add('allele_GRCh38:'+x['variant_id'],'gene:'+x['gene_id'],'AlphaGenome_website_gene_expression_score','Google DeepMind AlphaGenome Atlas','RNA-seq',x['biosample'],'','assay '+x['assay']+'; displayed raw '+x['displayed_raw_score']+'; website display '+x['displayed_other_score']+'; '+x['precision_and_limits'])
with open('outputs/stage5_preliminary_evidence_network.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(edges)
print('rows',len(edges),'new types',dict(collections.Counter(x['edge_type'] for x in edges if x['edge_type'].startswith('AlphaGenome') or x['edge_type']=='verified_candidate_allele_mapping')))
