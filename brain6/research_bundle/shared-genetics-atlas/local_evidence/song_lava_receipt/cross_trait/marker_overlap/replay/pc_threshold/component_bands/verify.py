import csv, hashlib, json, math
from pathlib import Path
p=Path(__file__).resolve().parent
d=json.loads((p/'projections.json').read_text()); prev=json.loads((p.parent/'results.json').read_text())
v=d['squared_projections']; assert len(v)==306 and all(x>=0 for x in v)
assert [r['K'] for r in prev]==d['threshold_K']
assert all(math.isclose(sum(v[:r['K']]),r['fi_delta_sumsq'],abs_tol=1e-12) for r in prev)
rows=[];start=0
for r in prev:
    end=r['K'];mean=sum(v[start:end])/(end-start)
    rows.append({'first_component':start+1,'last_component':end,'n':end-start,'mean_squared_projection':mean,
                 'sampling_term_at_upper_threshold':r['fi_sigma'],'mean_to_sampling_ratio':mean/r['fi_sigma']})
    start=end
with (p/'bands.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps({'bands':rows,'projection_sha256':hashlib.sha256((p/'projections.json').read_bytes()).hexdigest(),
                  'note':'Descriptive component sums and replay arithmetic; no statistical calibration.'},indent=2))
