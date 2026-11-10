import csv,collections,zipfile,math,pathlib,json,hashlib
with zipfile.ZipFile('work/stage4_hpa_cluster_metadata.zip') as z:
 meta=[x for x in csv.DictReader((s.decode() for s in z.open('rna_single_cell_clusters.tsv')),delimiter='\t') if x['Tissue']=='skeletal muscle']
chosen=[x for x in meta if x['Included in aggregation']=='yes'];assert len(meta)==28 and len(chosen)==25
ids=[x['Cluster'] for x in chosen];clusters={x['Cluster']:x for x in chosen};genes=collections.defaultdict(dict);sym={}
for x in csv.DictReader(open('work/stage4_hpa_skeletal_muscle_cluster.tsv'),delimiter='\t'):
 if x['Cluster'] not in clusters:continue
 g=x['Gene'];c=x['Cluster'];assert c not in genes[g];v=float(x['nCPM']);assert math.isfinite(v) and v>=0;genes[g][c]=math.log1p(v);sym[g]=x['Gene name']
assert all(len(x)==25 for x in genes.values()) and len(genes)==20162
p=pathlib.Path('work/stage4_magma/matrices/HPA_human_skeletal_muscle_GSE233882.txt')
with p.open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['GENE']+ids+['Average'])
 for g,v in sorted(genes.items()):
  a=[v[c] for c in ids];w.writerow([g]+a+[sum(a)/len(a)])
with open('outputs/stage4_hpa_skeletal_muscle_cluster_metadata.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=meta[0]);w.writeheader();w.writerows(meta)
candidate=collections.defaultdict(set)
for x in csv.DictReader(open('outputs/stage3_initial_variant_gene_table.csv')):
 g=x['gene ID'].split('.')[0]
 if g.startswith('ENSG'):candidate[g].add(x['Candidate gene'])
rows=[]
for g in sorted(candidate.keys() & genes.keys()):
 for c in ids:
  v=genes[g][c];rows.append({'gene_id':g,'hpa_gene_symbol':sym[g],'candidate_gene_symbols':'|'.join(sorted(candidate[g])),'cluster':c,'cell_type':clusters[c]['Cell type'],'cell_type_detail':clusters[c]['Cell type detail'],'cell_count':clusters[c]['Cell count'],'nCPM':math.expm1(v),'expression_evidence':'HPA aggregated cluster nCPM; not genetic association or causal mechanism'})
with open('outputs/stage4_hpa_muscle_candidate_expression.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
json.dump({'source_dataset':'GSE233882','source_expression_receipt':'stage4_hpa_muscle_source_receipt.json','metadata_url':'https://www.proteinatlas.org/download/tsv/rna_single_cell_clusters.tsv.zip','metadata_zip_sha256':hashlib.sha256(open('work/stage4_hpa_cluster_metadata.zip','rb').read()).hexdigest(),'source_gene_count':len(genes),'source_cluster_count':len(meta),'included_clusters':len(chosen),'excluded_clusters':[x['Cluster'] for x in meta if x['Included in aggregation']!='yes'],'transform':'log1p(nCPM); Average is equal mean across 25 included clusters, not a cell-weighted or measured tissue average','candidate_ids':len(candidate),'candidate_matched':len(candidate.keys() & genes.keys()),'candidate_missing':sorted(candidate.keys()-genes.keys()),'matrix_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},open('outputs/stage4_hpa_muscle_preparation_receipt.json','w'),indent=2)
print('genes',len(genes),'clusters',len(chosen),'candidate genes',len(candidate.keys() & genes.keys()),'candidate gene-cluster rows',len(rows))
