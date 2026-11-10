"""Variant-focused diagnostics on validation only; blind test stays unopened.

An activity-aware classifier can predict assay detectability mainly by recognizing
active regions. Report performance within experimentally active regions too.
The measured activity flag defines a reporting stratum, never an input feature.
"""
from pathlib import Path
import json,csv,numpy as np
from common import ROOT,WORK,OUT,STUDIES,metrics

rows=[]
for stem in ['salomon','siraj']:rows.extend(csv.DictReader((WORK/(stem+'_manifest_local.csv')).open()))
active=np.array([int(r['active_screen']) for r in rows],dtype='bool');result={}
files=[OUT/'stage1_validation_predictions.npz']
for representation in ['hyena_chromatin','motif','motif_chromatin']:
    prefix=representation+'_prior'
    if (OUT/(prefix+'_selection.json')).exists():files.append(OUT/(prefix+'_validation_predictions.npz'))
if (OUT/'full_avi_classical_finished.json').exists():files.append(OUT/'full_avi_classical_validation_predictions.npz')
for path in files:
    if not path.exists():continue
    with np.load(path) as z:
        ix=z['row_index'];y=z['label'];study=z['study'];names=[k for k in z.files if k not in ['row_index','label','study']]
        assert all(int(rows[i]['chromosome']) in [4,5] for i in ix)
        for name in names:
            pred=z[name];per={}
            for sid,label in enumerate(STUDIES):
                mask=(study==sid)&active[ix];allmask=study==sid;counts={'n':int(mask.sum()),'detected_calls':int(y[mask].sum()),'all_validation_detected_calls':int(y[allmask].sum()),'outside_active_detected_calls':int(y[allmask&~active[ix]].sum())}
                if mask.sum()>1 and len(np.unique(y[mask]))==2:counts['within_active_metrics']=metrics(y[mask],pred[mask])
                else:counts['within_active_metrics']=None
                per[label]=counts
            result[name]=per
receipt={'scope':'Validation chromosomes4-5 only; no test predictions read or model refitting','reason':'Assess variant-specific ranking within measured active regions, alongside all-variant detection performance','activity_flags':'Source Salomon significant reference/alternate region versus control; Siraj predefined source-row activity screen. Criteria differ and are not biological inactive/active truth.','activity_flag_used_as_feature':False,'models':result,'limitations':'Post-stage1 validation diagnostic specified before v3 test evaluation; conditional subsets are selected by laboratory measurements, not applicable to screening unmeasured variants without a predicted activity model.'}
(OUT/'validation_active_region_diagnostics.json').write_text(json.dumps(receipt,indent=2)+'\n');print('Within-active validation diagnostic models',len(result),flush=True)
