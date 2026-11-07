#!/usr/bin/env python3
import csv,gzip,hashlib,json,statistics
from pathlib import Path
root=Path(__file__).resolve().parent
old=Path('/Volumes/Extreme SSD/brain6-cfdr-mixed-signal-20261028')
new=Path('/Volumes/Brain6Work/cfdr-mixed-resume-20261028')
def table(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
rows=table(root/'runs.tsv')+table(root/'resume_runs.tsv')
assert len(rows)==36 and len({(r['seed'],r['rho'],r['method']) for r in rows})==36
hashes={};counts={}
for r in rows:
    markers,false,blocks,fb,shared,recovered=(int(r[k]) for k in ('markers','false_markers','blocks','false_blocks','shared_blocks','shared_blocks_detected'))
    assert 0<=false<=markers and 0<=fb<=blocks and 0<=recovered<=shared==22
    assert abs(float(r['fdp_markers'])-false/max(markers,1))<1e-12
    assert abs(float(r['fdp_blocks'])-fb/max(blocks,1))<1e-12
    key=(r['seed'],r['rho']);counts[key]=counts.get(key,0)+1
for (seed,rho),count in counts.items():
    assert count==4
    suffix=str(int(float(rho)*10));name=f'scores_{seed}_rho{suffix}.rds'
    candidates=[d/name for d in (old,new) if (d/name).exists()]
    assert len(candidates)==1
    path=candidates[0];h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    with gzip.open(path,'rb') as f:
        while f.read(1024*1024):pass
    hashes[str(path)]=h.hexdigest()
assert len(hashes)==9
summary=[]
for rho in ('0','0.3'):
 for method in ('smooth','unsmoothed','maxp_bonferroni','maxp_BY'):
    a=[r for r in rows if r['rho']==rho and r['method']==method]
    summary.append({'rho':rho,'method':method,'n':len(a),
      'mean_marker_fdp':statistics.mean(float(r['fdp_markers']) for r in a),
      'runs_false':sum(int(r['false_markers'])>0 for r in a),
      'mean_shared_blocks_detected':statistics.mean(int(r['shared_blocks_detected']) for r in a)})
receipt={'complete_datasets':9,'result_rows':36,'gzip_integrity':'pass','complete_file_sha256':hashes,
         'limit':'Checks table arithmetic and compressed-file integrity; does not decode RDS or independently verify resumed scores.',
         'incomplete_file':str(new/'scores_71005_rho3.rds'),'summary':summary}
(root/'resume_partial_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(summary,indent=2))
