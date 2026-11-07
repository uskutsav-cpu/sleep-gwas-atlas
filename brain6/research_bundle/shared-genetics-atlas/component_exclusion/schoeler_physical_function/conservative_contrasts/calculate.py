import itertools,json,math
from pathlib import Path
root=Path(__file__).resolve().parent
rows=json.loads((root.parent/'audit-20261029T161051-7498.json').read_text())['comparisons']
values={r['target']:r for r in rows}
output=[]
for a,b in itertools.combinations(('FI','FFS','HFRS'),2):
 d=values[a]['rg']-values[b]['rg']
 # LDSC logged four decimal places: allow worst-case rounding in two SEs and estimates.
 bound=values[a]['se']+values[b]['se']+0.0001
 z=max(0.,abs(d)-0.0001)/bound
 output.append({'contrast':a+' minus '+b,'difference':d,'worst_case_se_bound':bound,'conservative_z':z,'conservative_two_sided_p':math.erfc(z/math.sqrt(2)),'symmetric_approx_95_percent_interval':[d-0.0001-1.96*bound,d+0.0001+1.96*bound]})
(root/'summary.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
