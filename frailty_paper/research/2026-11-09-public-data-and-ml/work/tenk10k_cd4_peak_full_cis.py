exec(open('work/tenk10k_remote_parquet_probe.py').read().split("f=Remote(u,")[0])
import pyarrow.compute as pc,csv,json,hashlib
from pathlib import Path
f=Remote(u,3653870227);p=pq.ParquetFile(f);peak='chr5:136132117-136134284';groups=[];out=[]
for i in range(p.num_row_groups):
 s=p.metadata.row_group(i).column(0).statistics
 if s is None or not s.has_min_max or s.min<=peak<=s.max:groups.append(i)
for i in groups:
 t=p.read_row_group(i,columns=['phenotype_id']);mask=pc.equal(t['phenotype_id'],peak)
 if pc.any(mask).as_py():
  for r in p.read_row_group(i).filter(mask).to_pylist():r['row_group']=i;out.append(r)
assert out
path=Path('work/tenk10k_cd4_smad5_peak_full_cis.csv')
with path.open('w',newline='') as fp:
 w=csv.DictWriter(fp,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
receipt={'rows':len(out),'unique_variants':len(set(x['variant_id'] for x in out)),'row_groups_checked':groups,'object_size':3653870227,'source_repo_revision':'a8af2b48c3e9c3a9daac9fb8d5396b377fdb929b','file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes_read':f.total,'method':'Remote parquet footer row-group min/max and full requested peak records; full remote object checksum not verified'}
Path('outputs/stage6_tenk10k_cd4_peak_full_cis_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
