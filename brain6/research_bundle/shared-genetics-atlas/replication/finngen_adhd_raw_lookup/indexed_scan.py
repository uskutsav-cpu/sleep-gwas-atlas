#!/usr/bin/env python3
import csv, hashlib, json, math
from pathlib import Path
import pysam
root=Path(__file__).resolve().parent
atlas=root.parent.parent
url='https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_ADHD.gz?generation=1777989549884404'
idx=root/'source.tbi'; data=idx.read_bytes()
if hashlib.md5(data).hexdigest()!='dd3757119e402f1184174653496f1599':raise ValueError('index MD5 differs from pinned object')
with (atlas/'replication/six_lead_joint/six_leads.tsv').open() as f: leads=list(csv.DictReader(f,delimiter='\t'))
leads.append({'variant':'rs62264764','chr':'3','position38':'117920728','effect_allele':'A','other_allele':'G'})
fields=['variant','chr','position38','effect_allele','other_allele','status','source_ids','ref','alt','source_p','beta_aligned','se','af_alt']
out=[]
with pysam.TabixFile(url,index=str(idx)) as t:
 for v in leads:
  c=v['chr'];p=int(v['position38']); candidates=list(t.fetch(c,p-1,p)); matches=[]
  for record in candidates:
   x=record.split('\t')
   if len(x)!=13:raise ValueError('unexpected column count')
   if x[0]==c and int(x[1])==p and v['variant'] in x[4].split(',') and {x[2],x[3]}=={v['effect_allele'],v['other_allele']}:matches.append(x)
  row={k:v.get(k,'') for k in fields};row['status']='NOT_FOUND' if not matches else 'MULTIPLE_MATCHES' if len(matches)>1 else 'MATCHED'
  if len(matches)==1:
   x=matches[0];pv,beta,se,af=map(float,(x[6],x[8],x[9],x[10]))
   if not (0<=pv<=1 and math.isfinite(beta) and math.isfinite(se) and se>0 and 0<=af<=1):raise ValueError('invalid statistics')
   row.update(source_ids=x[4],ref=x[2],alt=x[3],source_p=x[6],beta_aligned=str(beta if v['effect_allele']==x[3] else -beta),se=x[9],af_alt=x[10])
  out.append(row)
with (root/'indexed_results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(root/'indexed_receipt.json').write_text(json.dumps({'source_url':url,'source_length_from_http_metadata':802195177,'source_md5_from_http_metadata':'e7aedc4602071840bbd5abe62d6fcdcc','source_full_download_verified':False,'index_generation':'1777989512207369','index_bytes':len(data),'index_md5':hashlib.md5(data).hexdigest(),'index_sha256':hashlib.sha256(data).hexdigest(),'pysam_version':pysam.__version__,'selection':'six frozen leads plus separately selected rs62264764','partial_random_access_only':True},indent=2)+'\n')
print((root/'indexed_results.tsv').read_text())
