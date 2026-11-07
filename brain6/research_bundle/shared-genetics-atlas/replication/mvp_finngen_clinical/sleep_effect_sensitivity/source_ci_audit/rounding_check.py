import itertools,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
receipt=json.loads((p.parents[2]/'mvp_six_leads/receipt.json').read_text())
raw=next(x for x in receipt['audit'] if x['variant']=='rs2582894')['mvp_exact_rows'][0]
vals=[float(raw[k]) for k in ('odds_ratio','ci_lower','ci_upper')];assert all(len(raw[k].split('.')[1])==4 for k in ('odds_ratio','ci_lower','ci_upper'))
pvals=[]
for signs in itertools.product((-1,1),repeat=3):
 odds,lo,hi=[v+s*0.00005 for v,s in zip(vals,signs)];se=(math.log(hi)-math.log(lo))/(2*1.95996398454);pvals.append(math.erfc(abs(math.log(odds)/se)/math.sqrt(2)))
out={'variant':'rs2582894','raw':{k:raw[k] for k in ('odds_ratio','ci_lower','ci_upper','p_value')},'rounding_assumption':'nearest 0.0001 for all three values','wald_95ci_p_min':min(pvals),'wald_95ci_p_max':max(pvals),'source_p':float(raw['p_value']),'note':'This assumes ordinary Wald 95% limits; the source procedure is not established.'}
(p/'rounding_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
