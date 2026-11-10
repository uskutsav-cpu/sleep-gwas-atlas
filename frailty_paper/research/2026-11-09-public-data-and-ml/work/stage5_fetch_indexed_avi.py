import csv, os, json, threading, concurrent.futures, time
from collections import defaultdict
import pysam
url='https://grr.iossifovlab.com/hg38/scores/AVI/alphagenome_variant_impact_score_snvs.tsv.gz'
index='/Volumes/Extreme SSD/codex-sleep-frailty-work/stage4_atlas_cache/avi_snvs.tsv.gz.tbi'
rows=list(csv.DictReader(open('outputs/stage5_score_queue.csv')))
groups=defaultdict(list)
for row in rows:
 if len(row['reference_allele'])==len(row['alternate_allele'])==1:
  groups[(row['chr'],int(row['position_GRCh38'])//100000)].append(row)
local=threading.local()
def fetch(item):
 (chrom,window), variants=item
 path=f'work/stage5_indexed_results/{chrom}_{window}.json'
 if os.path.exists(path): return path
 if not hasattr(local,'tabix'): local.tabix=pysam.TabixFile(url,index=index)
 wanted={(int(v['position_GRCh38']),v['reference_allele'],v['alternate_allele']):v['variant_key'] for v in variants}
 found={}
 positions=[p[0] for p in wanted]
 for line in local.tabix.fetch(chrom,min(positions)-1,max(positions)):
  fields=line.split('\t'); key=(int(fields[1]),fields[2],fields[3])
  if key in wanted: found[wanted[key]]={'raw_score':fields[4],'avi_phred':fields[5]}
 data={'chrom':chrom,'window':window,'requested':[v['variant_key'] for v in variants],'found':found}
 with open(path+'.tmp','w') as f: json.dump(data,f)
 os.replace(path+'.tmp',path)
 print(path,len(found),'/',len(variants),flush=True)
 return path
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for future in concurrent.futures.as_completed([pool.submit(fetch,x) for x in groups.items()]):
  try: future.result()
  except Exception as e: print('ERROR',repr(e),flush=True)
