exec(open('work/tenk10k_remote_parquet_probe.py').read().split("f=Remote(u,")[0])
f=Remote(u,3653870227);p=pq.ParquetFile(f);target='chr5:136132117-136134284';possible=[]
for i in range(p.num_row_groups):
 s=p.metadata.row_group(i).column(0).statistics
 if s is None or not s.has_min_max or s.min<=target<=s.max:possible.append(i)
print('possible',len(possible),possible,'bytes',f.total)
for i in range(p.num_row_groups):
 s=p.metadata.row_group(i).column(0).statistics
 if i%10==0:print(i,s.min,s.max)
