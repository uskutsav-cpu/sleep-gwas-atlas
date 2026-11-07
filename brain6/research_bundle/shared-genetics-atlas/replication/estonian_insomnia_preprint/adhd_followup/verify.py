#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[2];receipt=json.loads((p/'receipt.json').read_text());idx=atlas/'replication/finngen_adhd_raw_lookup/source.tbi'
with (p/'results.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==2 and {x['variant'] for x in rows}=={'rs1093463','rs1989903'}
positions={'rs1093463':('3',118025915,'C','T'),'rs1989903':('7',114497885,'A','G')}
with pysam.TabixFile(receipt['source_url'],index=str(idx)) as tab:
 for x in rows:
  ch,pos,ref,alt=positions[x['variant']];hits=[r.split('\t') for r in tab.fetch(ch,pos-1,pos)]
  hits=[r for r in hits if r[0]==ch and int(r[1])==pos and x['variant'] in r[4].split(',') and {r[2],r[3]}=={ref,alt}]
  assert len(hits)==1 and x['adhd_status']=='unique'
  q=hits[0];aligned=float(q[8])*(1 if q[3]==alt else -1)
  assert math.isclose(float(x['adhd_p']),float(q[6]),abs_tol=1e-12) and math.isclose(float(x['adhd_beta_alt']),aligned,abs_tol=1e-12)
  assert math.isclose(float(x['max_p']),max(float(x['sleep_p']),float(x['adhd_p'])),abs_tol=1e-12)
assert all(float(x['max_p'])>.05 for x in rows)
r={'variant_matches':len(rows),'both_max_p_above_005':True,'source_rows_refetched':True,'sha256':hashlib.sha256((p/'results.tsv').read_bytes()).hexdigest(),'limitations':'Checks indexed-source matching and arithmetic, not cohort independence, calibration after selection or a shared causal variant.'}
(p/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
