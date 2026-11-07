import csv, json, math, statistics
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'pairwise.tsv').open() as f: pairs=list(csv.DictReader(f,delimiter='\t'))
with (p/'marker_status.tsv').open() as f: markers=list(csv.DictReader(f,delimiter='\t'))
r=json.loads((p/'receipt.json').read_text());assert len(markers)==60 and sum(x['status']=='matched' for x in markers)==r['matched_markers']
assert len(pairs)==r['summaries']['all']['pairs']==1770
for x in pairs:
    a,b,d=map(float,(x['ukb_r'],x['eur_r'],x['abs_r2_diff']))
    # Do not silently clip approximate reference values outside [-1, 1].
    assert abs(a)<=1.01 and abs(b)<=1.000001 and math.isclose(d,abs(a*a-b*b),abs_tol=1e-12)
for name,subset in [('all',pairs),('within_100kb',[x for x in pairs if int(x['distance_bp'])<=100000]),('ukb_r2_ge_0.2',[x for x in pairs if float(x['ukb_r2'])>=.2])]:
    d=[float(x['abs_r2_diff']) for x in subset];s=r['summaries'][name]
    assert len(d)==s['pairs'] and math.isclose(statistics.median(d),s['median_abs_r2_diff'],abs_tol=1e-12)
    assert sum(x>.1 for x in d)==s['over_0.1']
print(json.dumps({'matched_markers':r['matched_markers'],'pairs':len(pairs),
 'ukb_abs_r_over_one':sum(abs(float(x['ukb_r']))>1 for x in pairs),'summaries':r['summaries'],
 'scope':'arithmetic and selected marker counts only; not an independent genotype or reference audit'},indent=2))
