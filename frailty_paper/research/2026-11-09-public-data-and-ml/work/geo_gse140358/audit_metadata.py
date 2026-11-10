import csv,gzip,collections,json,hashlib
from pathlib import Path
p=Path('work/geo_gse140358/matrix.txt.gz');out={'source':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140358/matrix/GSE140358_series_matrix.txt.gz','sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
with gzip.open(p,'rt') as f:
 for r in csv.reader(f,delimiter='\t'):
  if r and r[0]=='!Sample_geo_accession':out['sample_count']=len(r)-1
  if r and r[0]=='!Sample_characteristics_ch1' and r[1].startswith('condition:'):out['source_condition_labels']=dict(collections.Counter(x.split(': ',1)[1] for x in r[1:]))
  if r and r[0]=='!Series_platform_id':out['platform']=r[1]
  if r and r[0]=='!series_matrix_table_begin':break
assert out['sample_count']==sum(out['source_condition_labels'].values())
Path('outputs/stage6_gse140358_metadata_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
