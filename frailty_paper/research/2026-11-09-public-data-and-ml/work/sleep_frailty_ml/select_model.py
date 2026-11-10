"""Apply the pre-test validation-only model recommendation rule."""
from pathlib import Path
import json
OUT=Path(__file__).resolve().parents[2]/'outputs/sleep_frailty_ai_v1'
h=json.loads((OUT/'training_history.json').read_text());b=json.loads((OUT/'baseline_training_history.json').read_text())
scores={'kmer_logistic':max(r['validation_average_precision'] for r in b),'frozen_hyenadna':max(r['validation_average_precision'] for r in h if r['phase']=='frozen_head'),'finetuned_hyenadna':max(r['validation_average_precision'] for r in h if r['phase']=='finetune')}
scores['small_dna_cnn']=json.loads((OUT/'cnn_selection.json').read_text())['validation_average_precision']
selected=max(scores,key=scores.get)
result={'rule':'Largest validation average precision; test results not read by this selector','selected_model':selected,'validation_average_precision':scores,'scope':'Recommendation for selected neuronal reporter detected-call ranking only, not clinical sleep/frailty prediction'}
(OUT/'model_selection.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
