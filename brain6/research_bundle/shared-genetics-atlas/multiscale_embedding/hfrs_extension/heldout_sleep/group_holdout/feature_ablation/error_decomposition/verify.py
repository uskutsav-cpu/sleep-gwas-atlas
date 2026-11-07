#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
cells=read(p/'cell_deltas.tsv');groups=read(p/'group_deltas.tsv');r=json.loads((p/'receipt.json').read_text())
assert len(cells)==216 and len(groups)==48
for g in groups:
 xs=[x for x in cells if x['dropped_feature']==g['dropped_feature'] and x['held_group']==g['held_group']]
 assert len(xs)==int(g['n'])
 for x in xs:assert math.isclose(float(x['squared_error_dropped'])-float(x['squared_error_all_features']),float(x['difference_dropped_minus_all']),abs_tol=1e-12)
 assert math.isclose(sum(float(x['difference_dropped_minus_all']) for x in xs),float(g['sum_squared_error_difference']),abs_tol=1e-12)
 assert sum(float(x['difference_dropped_minus_all'])<0 for x in xs)==int(g['improved_cells'])
for feature in r['totals']:
 assert math.isclose(sum(float(g['sum_squared_error_difference']) for g in groups if g['dropped_feature']==feature),r['totals'][feature]['sum_delta'],abs_tol=1e-12)
subset=[x for x in cells if x['dropped_feature']=='parkinson' and x['held_group'] not in ('dashti_2019_duration','jones_accelerometer')];s=r['parkinson_posthoc_subgroups'];assert len(subset)==s['remaining_cells']==18
for key,col in [('remaining_rmse_all','squared_error_all_features'),('remaining_rmse_dropped','squared_error_dropped')]:assert math.isclose(math.sqrt(sum(float(x[col]) for x in subset)/len(subset)),s[key],abs_tol=1e-12)
special=sum(float(x['difference_dropped_minus_all']) for x in cells if x['dropped_feature']=='parkinson' and x['held_group'] in ('dashti_2019_duration','jones_accelerometer'))
assert math.isclose(special/r['totals']['parkinson']['sum_delta'],s['special_fraction_of_net_delta'],abs_tol=1e-12)
out={'cells_verified':len(cells),'groups_verified':len(groups),'group_totals_reconciled':True,'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('cell_deltas.tsv','group_deltas.tsv')},'limits':'Arithmetic check only, not independent validation, calibrated inference or causal attribution.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
