import csv,json,math,pathlib
p=pathlib.Path('outputs/shared-genetics-atlas/multiscale_embedding');a=list(csv.DictReader((p/'source_matrix.tsv').open(),delimiter='\t'));r=list(csv.DictReader((p/'reconstruction.tsv').open(),delimiter='\t'));v=list(csv.DictReader((p/'masked_predictions.tsv').open(),delimiter='\t'));s=json.loads((p/'summary.json').read_text());assert len(a)==len(r)==96 and len(v)==960 and len({(x['sleep_trait'],x['outcome']) for x in a})==96
for col in ('column_mean','rank2','additive'):
 val=math.sqrt(sum((float(x[col])-float(x['observed']))**2 for x in v)/len(v));k={'column_mean':'column_mean','rank2':'rank2_ALS','additive':'additive'}[col];assert abs(val-s['masked_cell_rmse'][k])<1e-12
assert sum(abs(float(x['rank2']))>1 for x in v)==4
print('PASS: 96 unique cells, 960 masked predictions, reported RMSE and out-of-range count')
