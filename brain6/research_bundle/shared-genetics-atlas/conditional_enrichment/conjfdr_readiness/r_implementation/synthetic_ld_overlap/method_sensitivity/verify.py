#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
def read(name):
 with (p/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
r=read('sensitivity.tsv');hits=read('default_hits.tsv')
assert {x['config'] for x in r}=={'default','no_adjust','no_smoothing','no_genomic_correction'}
assert all(x['status']=='OK' and int(x['finite'])==44000 and int(x['nonfinite'])==0 for x in r)
d={x['config']:x for x in r};assert int(d['default']['hits'])==1 and int(d['no_smoothing']['hits'])==0
assert len(hits)==1 and max(float(hits[0]['cfdr12']),float(hits[0]['cfdr21']))==float(hits[0]['conj_fdr'])
assert float(hits[0]['p_trait1'])>float(hits[0]['p_trait2'])
log=(p/'run.log').read_text();assert 'METHOD_SENSITIVITY_COMPLETED' in log
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
files=[p/x for x in ('PROTOCOL.md','run.R','run.log','sensitivity.tsv','default_hits.tsv')]
out={'sha256':{str(f):sha(f) for f in files},'summaries':r,'default_hit':hits[0],
 'default_hit_max_p_bonferroni':min(1,44000*max(float(hits[0]['p_trait1']),float(hits[0]['p_trait2']))),
 'warning_count_logged':log.count('Warning in ')}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
