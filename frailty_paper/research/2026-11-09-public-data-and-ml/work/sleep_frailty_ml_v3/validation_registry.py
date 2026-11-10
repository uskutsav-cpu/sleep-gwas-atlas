"""Discover completed default-fold candidates without touching the blind test."""
import json
from common import *

def completed_candidates():
    files=[]
    if (OUT/'stage1_selection.json').exists():files.append(OUT/'stage1_search_results.json')
    for p in OUT.glob('*_search_results.json'):
        if p.name=='stage1_search_results.json':continue
        prefix=p.name.removesuffix('_search_results.json')
        if (OUT/(prefix+'_selection.json')).is_file():files.append(p)
    if (OUT/'full_avi_classical_finished.json').is_file():files.append(OUT/'full_avi_classical_results.json')
    records={};prediction_files={}
    for path in sorted(set(files)):
        prefix=path.name.removesuffix('_search_results.json') if path.name.endswith('_search_results.json') else 'full_avi_classical'
        pred=OUT/(prefix+'_validation_predictions.npz')
        if not pred.is_file():continue
        with np.load(pred) as z:
            ix=z['row_index'];chrom=np.array([int(r['chromosome']) for stem in ['salomon','siraj'] for r in csv.DictReader((WORK/(stem+'_manifest_local.csv')).open())])
            if set(chrom[ix])!={4,5}:continue # Other folds are kept separate.
            names=set(z.files)-{'row_index','label','study','active'}
        rows=json.loads(path.read_text());assert {r['name'] for r in rows}<=names
        for row in rows:
            assert row['name'] not in records,row['name'];records[row['name']]=row;prediction_files[row['name']]=pred
    return records,prediction_files

if __name__=='__main__':
    r,p=completed_candidates();(OUT/'completed_validation_registry.json').write_text(json.dumps({'candidate_count':len(r),'models':list(r.values()),'scope':'Completed default-fold validation searches only; not final model selection or a blind test result'},indent=2)+'\n');print('Completed validation candidates',len(r))
