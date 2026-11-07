#!/usr/bin/env python3
import csv,hashlib,json,statistics
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'runs.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==20
assert {(int(r['seed']),float(r['rho'])) for r in rows}=={(s,rho) for s in range(69001,69011) for rho in (0,.3)}
assert all(int(r['n'])==44000 and int(r['nonnull_blocks'])==220 and int(r['pruning_iterations'])==5 for r in rows)
assert 'RECHECK_PASSED' in (p/'warning_recheck.log').read_text()
summary={}
for rho in (0,.3):
 q=[r for r in rows if float(r['rho'])==rho]; counts=[int(r['rejections']) for r in q]
 summary[str(rho)]={'runs':len(q),'runs_with_rejection':sum(x>0 for x in counts),'total_snp_rejections':sum(counts),
 'smallest_conjunction_value':min(float(r['min_conj']) for r in q),'median_run_minimum':statistics.median(float(r['min_conj']) for r in q)}
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as r:
  for c in iter(lambda:r.read(8*1024*1024),b''):h.update(c)
 return h.hexdigest()
root=Path('/Volumes/Brain6Work/cfdr-synthetic-ld-overlap-20261027/reference')
files=[p/'PROTOCOL.md',p/'run.R',p/'runs.tsv',p/'warning_recheck.R',p/'warning_recheck.log',root/'all_chr_perVariant.rds']+sorted(root.glob('chr*_LDpairs.rds'))
assert len(files)==28
result={'scenarios':summary,'files_sha256':{str(f):sha(f) for f in files},'notes':'No real data. All simulated conjunction hypotheses are null; count is per run and does not establish human-data FDR.'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary,indent=2))
