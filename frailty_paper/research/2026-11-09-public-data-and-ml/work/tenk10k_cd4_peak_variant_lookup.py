exec(open('work/tenk10k_remote_parquet_probe.py').read().split("f=Remote(u,")[0])
import pyarrow.compute as pc,json,hashlib
from pathlib import Path
f=Remote(u,3653870227);p=pq.ParquetFile(f);peak='chr5:136132117-136134284';variants=['5:136202653:A:C','5:136131620:C:G'];possible=[];found=[]
for i in range(p.num_row_groups):
 s=p.metadata.row_group(i).column(0).statistics
 if s is None or not s.has_min_max or s.min<=peak<=s.max:possible.append(i)
for i in possible:
 table=p.read_row_group(i,columns=['phenotype_id','variant_id']);mask=pc.and_(pc.equal(table['phenotype_id'],peak),pc.is_in(table['variant_id'],value_set=__import__('pyarrow').array(variants)))
 if pc.any(mask).as_py():
  full=p.read_row_group(i);matches=full.filter(mask).to_pylist()
  for row in matches:row['row_group']=i;found.append(row)
print('groups',possible,'downloaded_bytes',f.total,'matches',len(found));print(json.dumps(found,indent=2))
out={'repository_revision':'a8af2b48c3e9c3a9daac9fb8d5396b377fdb929b','source_file':'caQTL_common_TensorQTL/CD4_Naive/TenK10K.cis_qtl_pairs.chr5.parquet','repository_lfs_sha256':'379f55643c2958ca6c6ea7ad6b4488f4c266a93276406d2b24fea9baac8ed687','file_size':3653870227,'row_count':p.metadata.num_rows,'peak':peak,'requested_variants':variants,'candidate_row_groups':possible,'bytes_read_remote':f.total,'results':found,'caveat':'Indexed remote reads do not verify full-file hash; absence is conditional on metadata and exact representation.'}
Path('outputs/stage6_tenk10k_cd4_peak_variant_lookup.json').write_text(json.dumps(out,indent=2)+'\n')
