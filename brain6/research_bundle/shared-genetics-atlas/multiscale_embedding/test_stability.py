import csv,json,math,pathlib,collections
p=pathlib.Path('outputs/shared-genetics-atlas/multiscale_embedding')
def read(n):
 with (p/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
a=read('multistart_training_objectives.tsv');b=read('multistart_predictions.tsv');s=json.loads((p/'multistart_summary.json').read_text());groups=collections.defaultdict(list)
for x in a:groups[(x['repeat'],x['fold'])].append(x)
assert len(groups)==50 and len(a)==250 and len(b)==960
for g in groups.values():
 chosen=[x for x in g if x['selected']=='True'];assert len(chosen)==1 and float(chosen[0]['training_objective'])==min(float(x['training_objective']) for x in g)
for k in ('column_mean','original_rank2','additive','multistart_rank2'):
 v=math.sqrt(sum((float(r[k])-float(r['observed']))**2 for r in b)/len(b));assert abs(v-s['rmse'][k])<1e-12
print('PASS: 50 split-specific training-only selections, 250 fits, and 960 evaluated cells')
