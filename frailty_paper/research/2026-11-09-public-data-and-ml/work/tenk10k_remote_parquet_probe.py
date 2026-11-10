import io,requests,pyarrow.parquet as pq
class Remote(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.s=requests.Session();self.total=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos)
  if n<=0:return b''
  r=self.s.get(self.url,headers={'Range':f'bytes={self.pos}-{self.pos+n-1}'},timeout=60);r.raise_for_status()
  if r.status_code!=206 or len(r.content)!=n:raise IOError(f'invalid range: status {r.status_code}, length {len(r.content)}, expected {n}')
  self.pos+=n;self.total+=n;return r.content
u='https://huggingface.co/datasets/anglixue/TenK10K_multiome/resolve/a8af2b48c3e9c3a9daac9fb8d5396b377fdb929b/caQTL_common_TensorQTL/CD4_Naive/TenK10K.cis_qtl_pairs.chr5.parquet'
f=Remote(u,3653870227);p=pq.ParquetFile(f);print('rows',p.metadata.num_rows,'rowgroups',p.num_row_groups,'schema',p.schema.names,'bytes_fetched',f.total)
for i in range(min(3,p.num_row_groups)):
 g=p.metadata.row_group(i)
 print('rowgroup',i,'rows',g.num_rows,[(g.column(j).path_in_schema,str(g.column(j).statistics)[:300]) for j in range(min(2,g.num_columns))])
