"""Fixed ensemble inference and outcome-free confidence input construction."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from matched_coverage_head import MatchedCoverageHead,torch
from safetensors.torch import load_file
DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation'

def aggregate(outputs):
    keys=['magnitude','signed_effect','activity_reference','activity_alternate'];result={k:np.mean([p[k] for p in outputs],axis=0) for k in keys};result['ensemble_sd']=np.std([p['magnitude'] for p in outputs],axis=0)
    if 'observation_index' in outputs[0]:
        ids=outputs[0]['observation_index'];assert all(np.array_equal(p['observation_index'],ids) for p in outputs);result['observation_index']=ids
    return result

def load_ensemble(study,arm,other_calibration=False):
    settings=json.loads((FINAL/'final_settings_lock.json').read_text());records=[];models=[];predictions=[]
    for seed in settings['seeds']:
        r=json.loads((FINAL/'fits'/f'fixed_{study}_{arm}_{seed}_record.json').read_text());assert r['settings_lock_sha256']==sha(FINAL/'final_settings_lock.json');assert sha(r['checkpoint'])==r['checkpoint_sha256']
        key='other_calibration_predictions' if other_calibration else 'predictions';assert sha(r[key])==r[key+'_sha256']
        model=MatchedCoverageHead(r['cell_medians'],r['line_medians'],r['line_parent'],seed,arm=='AVI');model.load_state_dict(load_file(r['checkpoint']));model.eval();records.append(r);models.append(model)
        with np.load(r[key]) as z:predictions.append({k:z[k] for k in z.files})
    for r in records[1:]:
        for key in ['cell_medians','line_medians','supported_cells','supported_lines','training_cell_source_counts','training_line_source_counts','feature_mean','feature_sd','AVI_transform']:assert r[key]==records[0][key]
    return records[0],models,aggregate(predictions)

def confidence_inputs(pred,x,cell,line,metadata,exact_line_known=None):
    cell=np.asarray(cell,dtype='int64');line=np.asarray(line,dtype='int64');known=np.ones(len(cell),dtype=bool) if exact_line_known is None else np.asarray(exact_line_known,dtype=bool)
    cm=np.asarray(metadata['cell_medians'])[cell];lm=np.asarray(metadata['line_medians'])[line];cs=np.asarray(metadata['supported_cells'])[cell];ls=np.asarray(metadata['supported_lines'])[line]&known
    baseline=np.where(ls,lm,cm);counts=np.asarray(metadata['training_cell_source_counts'])[cell];lcounts=np.asarray(metadata['training_line_source_counts'])[line]*known
    features=np.column_stack([pred['magnitude'],pred['ensemble_sd'],abs(pred['signed_effect']),(pred['activity_reference']+pred['activity_alternate'])*.5,cm,baseline,counts,lcounts,cs,ls,x,np.eye(len(metadata['cell_medians']),dtype='float32')[cell]])
    assert np.isfinite(features).all();return features.astype('float32'),baseline,cm,cs,ls,counts

def predict_inputs(models,metadata,bank,x,cell,line,known_line,avi=None):
    nc=len(metadata['cell_medians']);nl=len(metadata['line_medians']);cs=np.asarray(metadata['supported_cells'])[cell];ls=np.asarray(metadata['supported_lines'])[line]&known_line
    cf=np.eye(nc,dtype='float32')[cell]*cs[:,None];lf=np.eye(nl,dtype='float32')[line]*ls[:,None];normalized=np.clip((x-np.asarray(metadata['feature_mean']))/np.asarray(metadata['feature_sd']),-12,12).astype('float32')
    av=None
    if metadata['arm']=='AVI':
        assert avi is not None;t=metadata['AVI_transform'];av=avi.copy();av[:,:2]=(np.where(np.isfinite(av[:,:2]),av[:,:2],t['median'])-t['mean'])/t['sd'];av=np.clip(av,-12,12).astype('float32')
    results=[]
    for model in models:
        output={k:np.empty(len(cell),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for i in range(0,len(cell),2048):
                ix=slice(i,i+2048);got=model(torch.from_numpy(normalized[ix]),torch.from_numpy(bank[ix,:,0].copy()),torch.from_numpy(bank[ix,:,1].copy()),torch.from_numpy(cell[ix].astype('int64')),torch.from_numpy(cf[ix]),torch.from_numpy(lf[ix]),torch.from_numpy(av[ix]) if av is not None else None)
                for k in output:output[k][ix]=got[k].numpy()
        results.append(output)
    return aggregate(results)

def weighted_quantile(values,groups,q=.9):
    values=np.asarray(values);groups=np.asarray(groups);assert len(values) and np.isfinite(values).all();_,inverse,counts=np.unique(groups,return_inverse=True,return_counts=True);weights=1./counts[inverse];order=np.argsort(values);j=np.searchsorted(np.cumsum(weights[order]),q*weights.sum(),side='left');return float(values[order[min(j,len(order)-1)]])
