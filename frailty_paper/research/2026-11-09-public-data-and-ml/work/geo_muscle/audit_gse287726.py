import csv,gzip,json,hashlib
from pathlib import Path
m=Path('work/geo_muscle/GSE287726.txt.gz');c=Path('work/geo_muscle/GSE287726_counts.tsv.gz');metadata=[]
with gzip.open(m,'rt') as f:
 for row in csv.reader(f,delimiter='\t'):
  if row and row[0]=='!Sample_title':titles=row[1:]
  if row and row[0]=='!Sample_characteristics_ch1':metadata.append(sorted(set(x.split(': ',1)[0] for x in row[1:] if ': ' in x)))
  if row and row[0]=='!series_matrix_table_begin':break
genes={'FOXP2':'ENSG00000128573','SMAD5':'ENSG00000113658','RBM6':'ENSG00000004534','MEIS1':'ENSG00000143995'};found={}
with gzip.open(c,'rt') as f:
 r=csv.reader(f,delimiter='\t');head=next(r)[1:];assert len(head)==70 and len(set(head))==70 and set(head)==set(titles)
 for row in r:
  if not row:continue
  ident=row[0].split('|')[0].split('.')[0]
  for g,i in genes.items():
   if ident==i:found[g]=row[0]
out={'n':len(titles),'characteristics_keys':metadata,'candidate_gene_rows':found,'metadata_sha256':hashlib.sha256(m.read_bytes()).hexdigest(),'counts_sha256':hashlib.sha256(c.read_bytes()).hexdigest(),'limitations':'Public GEO sample characteristics only identify tissue and age 70+; frailty per sample and sex not established. No candidate group comparison performed.'}
Path('outputs/stage6_gse287726_data_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
