"""Recompute frozen evaluation metrics; no model fitting or new inference."""
from common import *
def main():
    report=json.loads((OUT/'reserved_evaluation_report.json').read_text());assert report['prediction_sha256']==sha(OUT/'reserved_test_predictions.npz');exposure=json.loads((OUT/'reserved_exposure_audit.json').read_text());manifest=[]
    for stem in ['salomon','siraj']:
        p=WORK/(stem+'_manifest_local.csv');assert exposure['source_manifest_sha256'][stem]==sha(p);manifest.extend(csv.DictReader(p.open()))
    y=np.array([int(r['label']) for r in manifest]);study=np.array([STUDIES.index(r['study']) for r in manifest]);active=np.array([int(r['active_screen']) for r in manifest]);checks=[]
    with np.load(OUT/'reserved_test_predictions.npz') as z:
        ix=z['row_index'];strict=np.isin(ix,exposure['retained_Siraj_global_row_indices'])
        for s,name in enumerate(STUDIES):
            for stratum,suffix in [('all_calls','__overall'),('measured_active_regions','__active_variant')]:
                use=(study[ix]==s)&(strict if s==1 else True)&((active[ix]==1) if stratum!='all_calls' else True)
                for model,expected in report['results'][name]['strata'][stratum]['point_metrics'].items():
                    actual=metrics(y[ix[use]],z[model+suffix][use]);errors={k:abs(actual[k]-expected[k]) for k in ['auroc','average_precision','prevalence','brier_uncalibrated']};assert max(errors.values())<1e-12 and actual['n']==expected['n'];checks.append({'study':name,'stratum':stratum,'model':model,'errors':errors})
    (OUT/'frozen_metric_reproduction_check.json').write_text(json.dumps({'checks':checks,'prediction_sha256':report['prediction_sha256'],'purpose':'Frozen-prediction metric reconstruction; no new model inference, evaluation selection, fitting or new independent validation.'},indent=2)+'\n');print('Reproduced frozen metric rows',len(checks),flush=True)
if __name__=='__main__':main()
