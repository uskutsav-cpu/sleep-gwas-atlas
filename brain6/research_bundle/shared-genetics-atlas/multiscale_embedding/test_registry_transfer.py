import pathlib,csv,collections,math,json,hashlib
p=pathlib.Path('outputs/shared-genetics-atlas/multiscale_embedding')
with (p/'registry_transfer_predictions.tsv').open() as f:a=list(csv.DictReader(f,delimiter='\t'))
with (p/'registry_transfer_summary.tsv').open() as f:b=list(csv.DictReader(f,delimiter='\t'))
assert len(a)==80 and len(b)==8
for x in b:
 y=[r for r in a if r['outcome']==x['outcome'] and r['method']==x['method']];assert len(y)==10 and len({r['repeat'] for r in y})==10
 for suffix in ('UKB','R13'):
  actual=math.sqrt(sum(float(r['error_vs_'+suffix])**2 for r in y)/len(y));assert abs(actual-float(x['RMSE_vs_'+suffix]))<1e-12
print('PASS: two targets, four methods, ten repeated masks and all RMSE summaries')
