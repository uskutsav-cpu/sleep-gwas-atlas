"""Fixed inference and uncertainty features, with no fresh outcomes."""
from v9_common import *
from feature_banks import global_bank,local_bank
from line_aware_network import LineAwareEffectNetwork,torch
from safetensors.torch import load_file

FINAL=ROOT/'outputs/sleep_frailty_ai_v9_validation'
STATIC_NAMES=['global_reference_mean','global_alternate_mean','global_reference_sd','global_alternate_sd',
              'global_absolute_difference_mean','global_absolute_difference_max',
              'local_absolute_difference_mean','local_absolute_difference_max',
              'local_reference_mean','local_alternate_mean','local_reference_sd','local_alternate_sd']


def static_features(g,l):
    return np.column_stack([g[:,0].mean(1),g[:,1].mean(1),g[:,0].std(1),g[:,1].std(1),
                            abs(g[:,1]-g[:,0]).mean(1),abs(g[:,1]-g[:,0]).max(1),
                            abs(l[:,1]-l[:,0]).mean(1),abs(l[:,1]-l[:,0]).max(1),
                            l[:,0].mean(1),l[:,1].mean(1),l[:,0].std(1),l[:,1].std(1)]).astype('float32')


def load_ensemble(study,arm):
    lock=json.loads((FINAL/'final_settings_lock.json').read_text());records=[];models=[];predictions=[]
    for seed in lock['seeds']:
        path=FINAL/'fits'/f'fixed_{study}_{arm}_{seed}_record.json';r=json.loads(path.read_text())
        assert r['settings_lock_sha256']==sha(FINAL/'final_settings_lock.json')
        assert sha(r['checkpoint'])==r['checkpoint_sha256'] and sha(r['predictions'])==r['predictions_sha256']
        model=LineAwareEffectNetwork(r['cell_medians'],r['line_medians'],r['line_parent_cell'],r['recipe'],seed,arm=='AVI')
        model.load_state_dict(load_file(r['checkpoint']));model.eval();records.append(r);models.append(model)
        with np.load(r['predictions']) as z:predictions.append({k:z[k] for k in z.files})
    for r in records[1:]:
        for key in ['cell_medians','line_medians','supported_cells','supported_lines','training_cell_source_counts','training_line_source_counts','AVI_transform']:
            assert r[key]==records[0][key],key
    for p in predictions[1:]:assert np.array_equal(p['observation_index'],predictions[0]['observation_index'])
    return records[0],models,aggregate(predictions)


def aggregate(predictions):
    out={k:np.mean([p[k] for p in predictions],axis=0) for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
    out['ensemble_sd']=np.std([p['magnitude'] for p in predictions],axis=0)
    if 'observation_index' in predictions[0]:out['observation_index']=predictions[0]['observation_index']
    return out


def confidence_features(pred,static,cell,line,metadata):
    nc=len(metadata['cell_medians']);cmed=np.asarray(metadata['cell_medians'])[cell];lmed=np.asarray(metadata['line_medians'])[line]
    cs=np.asarray(metadata['supported_cells'])[cell];ls=np.asarray(metadata['supported_lines'])[line]
    baseline=np.where(ls,lmed,cmed)
    counts=np.asarray(metadata['training_cell_source_counts'])[cell];line_counts=np.asarray(metadata['training_line_source_counts'])[line]
    x=np.column_stack([pred['magnitude'],pred['ensemble_sd'],abs(pred['signed_effect']),
                       (pred['activity_reference']+pred['activity_alternate'])*.5,cmed,baseline,counts,line_counts,cs,ls,static,
                       np.eye(nc,dtype='float32')[cell]])
    names=['predicted_magnitude','ensemble_sd','absolute_predicted_allele_activity_difference','predicted_mean_allele_activity',
           'training_cell_median','training_line_or_cell_median','same_cell_training_studies','same_line_training_studies','known_cell','known_line']+STATIC_NAMES+['cell_'+s for s in metadata['cell_names']]
    assert np.isfinite(x).all()
    return x.astype('float32'),baseline,cmed,cs,ls,counts,names


def predict_blind(models,metadata,features,cell,line,avi):
    results=[];nc=len(metadata['cell_medians']);nl=len(metadata['line_medians'])
    cf=np.eye(nc,dtype='float32')[cell]*np.asarray(metadata['supported_cells'],dtype='float32')[cell,None]
    lf=np.eye(nl,dtype='float32')[line]*np.asarray(metadata['supported_lines'],dtype='float32')[line,None]
    av=None
    if metadata['arm']=='AVI':
        t=metadata['AVI_transform'];a=np.array(avi,copy=True);a[:,:2]=(np.where(np.isfinite(a[:,:2]),a[:,:2],t['median'])-t['mean'])/t['sd'];av=np.clip(a,-12,12).astype('float32')
    for model in models:
        out={k:np.empty(len(cell),dtype='float32') for k in ['magnitude','signed_effect','activity_reference','activity_alternate']}
        with torch.no_grad():
            for start in range(0,len(cell),2048):
                ix=slice(start,start+2048);h=torch.from_numpy(features[ix]);got=model(h[:,0,0],h[:,1,0],h[:,0,1],h[:,1,1],torch.from_numpy(cell[ix].astype('int64')),torch.from_numpy(line[ix].astype('int64')),torch.from_numpy(cf[ix]),torch.from_numpy(lf[ix]),torch.from_numpy(av[ix]) if av is not None else None)
                for k in out:out[k][ix]=got[k].numpy()
        results.append(out)
    return aggregate(results)


def weighted_quantile(values,groups,q=.9):
    unique,inverse,counts=np.unique(groups,return_inverse=True,return_counts=True);weights=1./counts[inverse]
    order=np.argsort(values);i=np.searchsorted(np.cumsum(weights[order]),q*weights.sum(),side='left');return float(values[order[min(i,len(order)-1)]])


def family(cell_name):
    if cell_name in ['neural_progenitor','induced_excitatory_neuron','neuronal_tumor','glioma']:return 'neural'
    if cell_name=='hematopoietic':return 'blood'
    return 'other'


def calibration_factors(ratios,cell,groups,cell_names):
    assert len(ratios) and np.isfinite(ratios).all()
    global_factor=weighted_quantile(ratios,groups);factors={};details={}
    fam=np.array([family(cell_names[c]) for c in cell])
    for c,name in enumerate(cell_names):
        ix=cell==c;scope='cell'
        if ix.sum()<100 or len(np.unique(groups[ix]))<20:
            ix=fam==family(name);scope='family'
            if ix.sum()<200 or len(np.unique(groups[ix]))<30:ix=np.ones(len(ratios),dtype=bool);scope='global'
        factors[name]=max(weighted_quantile(ratios[ix],groups[ix]),1e-4)
        details[name]={'factor':factors[name],'scope':scope,'observations':int(ix.sum()),'regions':len(np.unique(groups[ix]))}
    return {'target_coverage':.9,'global_factor':global_factor,'by_cell':factors,'details':details,
            'claim':'Empirical region-balanced source-omission calibration; no universal cross-study coverage guarantee.'}
