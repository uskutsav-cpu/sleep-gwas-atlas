#!/usr/bin/env python3
"""Descriptive bulk read-count detection; no association inference."""
import csv,gzip,json,hashlib,collections
from pathlib import Path
M=Path('work/GSE167186_GPL20301_series_matrix.txt.gz'); C=Path('work/GSE167186_counts.csv.gz'); selected=['FOXP2','SMAD5','RBM6','MEIS1']
records={}
with gzip.open(M,'rt') as f:
 for line in f:
  if line.startswith('!Sample_title\t'): records['sample'] = next(csv.reader([line],delimiter='\t'))[1:]
  if line.startswith('!Sample_characteristics_ch1\t'):
   vals=next(csv.reader([line],delimiter='\t'))[1:]
   if all(v.startswith('group: ') for v in vals): records['group']=[v.removeprefix('group: ') for v in vals]
assert len(records['sample'])==len(records['group'])==76 and len(set(records['sample']))==76
with gzip.open(C,'rt',newline='') as f:
 reader=csv.reader(f);head=next(reader);assert head[1:]==records['sample'];found={g:[] for g in selected}; total=0
 for row in reader:
  total+=1
  if row[0] in found:found[row[0]].append([int(x) for x in row[1:]])
result=[]
for gene in selected:
 matches=found[gene]
 for group in sorted(set(records['group'])):
  idx=[i for i,g in enumerate(records['group']) if g==group]
  values=[matches[0][i] for i in idx] if len(matches)==1 else []
  result.append({'gene':gene,'group':group,'n_samples':len(idx),'matching_gene_rows':len(matches),'samples_count_at_least_1':sum(x>=1 for x in values) if values else '', 'samples_count_at_least_10':sum(x>=10 for x in values) if values else ''})
with Path('outputs/stage6_gse167186_bulk_candidate_detection.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=result[0].keys());w.writeheader();w.writerows(result)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
receipt={'group_counts':dict(collections.Counter(records['group'])),'gene_rows':{x:len(found[x]) for x in selected},'genes_scanned':total,'input_sha256':{'matrix':digest(M),'counts':digest(C)},'raw_read_counts_not_normalized':True,'interpretation':'Descriptive bulk measurement detection only; no cell type, differential expression, genetic association, clinical validation or causal inference.'}
Path('outputs/stage6_gse167186_bulk_candidate_detection.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
