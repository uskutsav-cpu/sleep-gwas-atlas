from pathlib import Path
import csv,json
import numpy as np
OUT=Path(__file__).resolve().parents[2]/'outputs/sleep_frailty_ai_v1'
result=json.loads((OUT/'test_results.json').read_text());result['models']['small_dna_cnn']=json.loads((OUT/'cnn_test_results.json').read_text());result['cnn_selection']=json.loads((OUT/'cnn_selection.json').read_text());(OUT/'test_results.json').write_text(json.dumps(result,indent=2)+'\n')
data=np.load(OUT/'cnn_test_predictions.npz');lookup={int(i):float(p) for i,p in zip(data['row_index'],data['predictions'])}
with (OUT/'test_predictions.csv').open() as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
assert len(rows)==len(lookup)==3542
if 'small_dna_cnn' not in fields:fields.append('small_dna_cnn')
for row in rows:row['small_dna_cnn']=lookup[int(row['row_index'])]
with (OUT/'test_predictions.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
print('Merged CNN comparison on identical final-test row identifiers')
