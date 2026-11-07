#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[2];previous=atlas/'replication/finngen_adhd_raw_lookup';idx=previous/'source.tbi'
assert hashlib.md5(idx.read_bytes()).hexdigest()=='dd3757119e402f1184174653496f1599'
source=json.loads((previous/'indexed_receipt.json').read_text());url=source['source_url']
with (p.parent/'matched_candidates.tsv').open() as f: variants=list(csv.DictReader(f,delimiter='\t'))
assert {x['rsid'] for x in variants}=={'rs1093463','rs1989903'}
out=[]
with pysam.TabixFile(url,index=str(idx)) as tab:
 for v in variants:
  ch=v['chromosome'];pos=int(v['position38']);candidates=[]
  for line in tab.fetch(ch,pos-1,pos):
   x=line.split('\t');assert len(x)==13
   if x[0]==ch and int(x[1])==pos and v['rsid'] in x[4].split(',') and {x[2],x[3]}=={v['allele0'],v['allele1']}:candidates.append(x)
  result={'variant':v['rsid'],'position38':pos,'alt':v['allele1'],'sleep_p':float(v['p']),'adhd_status':'absent' if not candidates else 'unique' if len(candidates)==1 else 'ambiguous','adhd_p':'','adhd_beta_alt':'','adhd_se':'','max_p':''}
  if len(candidates)==1:
   x=candidates[0];pv,b,se,af=map(float,(x[6],x[8],x[9],x[10]));assert 0<=pv<=1 and math.isfinite(b) and se>0 and 0<=af<=1
   result.update(adhd_p=pv,adhd_beta_alt=b if x[3]==v['allele1'] else -b,adhd_se=se,max_p=max(float(v['p']),pv))
  out.append(result)
with (p/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'source_url':url,'pinned_index_md5':hashlib.md5(idx.read_bytes()).hexdigest(),'pinned_index_sha256':hashlib.sha256(idx.read_bytes()).hexdigest(),'previous_full_source_audit_sha256':hashlib.sha256((previous/'whole_source_hash_receipt.json').read_bytes()).hexdigest(),'selection':'two post-selected variants appearing in Estonian insomnia supplement','results':out,'limitations':'Indexed lookup, post-selected phenotype and variants. No proof of cohort or family independence, causal colocalization, or replicated sleep–ADHD mechanisms.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
