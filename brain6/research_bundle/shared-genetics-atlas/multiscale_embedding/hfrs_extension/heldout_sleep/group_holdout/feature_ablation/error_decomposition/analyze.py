#!/usr/bin/env python3
import collections,csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;source=p.parent/'predictions.tsv';summary=p.parent/'summary.tsv'
with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with summary.open() as f:previous=list(csv.DictReader(f,delimiter='\t'))
index={(x['dropped_feature'],x['sleep_trait'],x['target']):x for x in rows};assert len(index)==252
features=['adhd','mdd','scz','bipolar','alz','parkinson'];baseline=[x for x in rows if x['dropped_feature']=='none'];assert len(baseline)==36
out=[]
for dropped in features:
 for a in baseline:
  b=index[(dropped,a['sleep_trait'],a['target'])];assert b['held_group']==a['held_group'] and float(b['observed'])==float(a['observed'])
  y=float(a['observed']);e0=(float(a['prediction'])-y)**2;e1=(float(b['prediction'])-y)**2
  out.append({'dropped_feature':dropped,'sleep_trait':a['sleep_trait'],'held_group':a['held_group'],'target':a['target'],'squared_error_all_features':e0,'squared_error_dropped':e1,'difference_dropped_minus_all':e1-e0})
groups=collections.defaultdict(list)
for x in out:groups[(x['dropped_feature'],x['held_group'])].append(x)
g=[]
for (feature,group),rs in sorted(groups.items()):
 g.append({'dropped_feature':feature,'held_group':group,'n':len(rs),'sum_squared_error_difference':sum(x['difference_dropped_minus_all'] for x in rs),'improved_cells':sum(x['difference_dropped_minus_all']<0 for x in rs),'worse_cells':sum(x['difference_dropped_minus_all']>0 for x in rs)})
totals={}
for feature in features:
 items=[x for x in out if x['dropped_feature']==feature];assert len(items)==36
 prior0=next(x for x in previous if x['target']=='ALL' and x['dropped_feature']=='none');prior1=next(x for x in previous if x['target']=='ALL' and x['dropped_feature']==feature)
 delta=sum(x['difference_dropped_minus_all'] for x in items);assert math.isclose(delta/36,float(prior1['rmse'])**2-float(prior0['rmse'])**2,abs_tol=1e-12)
 totals[feature]={'sum_delta':delta,'improved_cells':sum(x['difference_dropped_minus_all']<0 for x in items),'worse_cells':sum(x['difference_dropped_minus_all']>0 for x in items)}
def write(path,rows):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write(p/'cell_deltas.tsv',out);write(p/'group_deltas.tsv',g)
rest=[x for x in out if x['dropped_feature']=='parkinson' and x['held_group'] not in ('dashti_2019_duration','jones_accelerometer')]
assert len(rest)==18
special=sum(x['difference_dropped_minus_all'] for x in out if x['dropped_feature']=='parkinson' and x['held_group'] in ('dashti_2019_duration','jones_accelerometer'))
subgroup={'special_group_sum_delta':special,'special_fraction_of_net_delta':special/totals['parkinson']['sum_delta'],'remaining_cells':len(rest),'remaining_rmse_all':math.sqrt(sum(x['squared_error_all_features'] for x in rest)/len(rest)),'remaining_rmse_dropped':math.sqrt(sum(x['squared_error_dropped'] for x in rest)/len(rest))}
receipt={'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'cell_rows':len(out),'group_rows':len(g),'totals':totals,'parkinson_posthoc_subgroups':subgroup,'limitations':'Post-hoc deterministic decomposition of dependent held-out predictions, not external validation, inference, causal feature importance or model selection.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
