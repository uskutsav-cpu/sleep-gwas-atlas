"""Independent raw-effect MAE replay from saved predictions and immutable data."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ml_resumable'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recipe',default='frozen_reader');parser.add_argument('--source',default='none');parser.add_argument('--epoch',type=int,default=1);args=parser.parse_args()
    directory=OUT/'development';stem=f'{args.recipe}_{args.source}_42'
    meta=json.loads((directory/(stem+'_started.json')).read_text())
    receipt=json.loads((directory/f'{stem}_epoch{args.epoch}_receipt.json').read_text())
    predpath=directory/f'{stem}_epoch{args.epoch}_predictions.npz';assert sha(predpath)==receipt['prediction_sha256']
    checkpoint=directory/f'{stem}_epoch{args.epoch}.safetensors';assert sha(checkpoint)==receipt['checkpoint_sha256']
    plan=json.loads((OUT/'finite_stable_native_plan.json').read_text());data=ROOT/'outputs/sleep_frailty_ml_next/coverage_training_arrays.npz'
    assert sha(data)==plan['data_sha256']
    with np.load(data,allow_pickle=False) as arrays:
        y=arrays['effect'];cell=arrays['cell'];line=arrays['line'];split=arrays['split'];study=arrays['study'];group=arrays['group']
        names=arrays['cell_names'].tolist();study_names=arrays['study_names'].tolist()
    train=split==0
    if args.source=='none':test=split==1
    else:
        sid=study_names.index(args.source);excluded=np.unique(group[study==sid]);train&=(study!=sid)&~np.isin(group,excluded);test=(study==sid)&(split!=2)
    tr,te=np.flatnonzero(train),np.flatnonzero(test)
    assert not np.intersect1d(group[tr],group[te]).size
    assert hashlib.sha256(tr.tobytes()).hexdigest()==meta['training_index_sha256']
    assert hashlib.sha256(te.tobytes()).hexdigest()==meta['prediction_index_sha256']
    with np.load(predpath,allow_pickle=False) as values:
        assert np.array_equal(values['observation_index'],te);prediction=values['magnitude'].astype('float64')
    assert np.isfinite(prediction).all() and np.all(prediction>=0) and np.all(y[te]>=0)
    cm=np.array(meta['cell_medians']);lm=np.array(meta['line_medians']);supported=np.array(meta['supported_cells']);lines=np.array(meta['supported_lines'])
    baseline=np.where(lines[line[te]],lm[line[te]],cm[cell[te]])
    curve=json.loads((directory/(stem+'_curve.json')).read_text())[args.epoch-1]
    cells={};skills=[];errors=[]
    for c in np.unique(cell[te]):
        ix=cell[te]==c;actual=y[te][ix].astype('float64')
        mae=float(np.mean(np.abs(actual-prediction[ix])));be=float(np.mean(np.abs(actual-baseline[ix])));skill=1-mae/be
        reported=curve['cells'][names[c]]
        errors.extend([abs(mae-reported['MAE']),abs(be-reported['line_baseline_MAE']),abs(skill-reported['MAE_skill_vs_line_baseline'])])
        cells[names[c]]={'rows':int(ix.sum()),'raw_observed_effect_MAE':mae,'training_typical_effect_MAE':be,'skill':skill,'supported':bool(supported[c])}
        if supported[c]:skills.append(skill)
    macro=float(np.mean(skills));errors.append(abs(macro-curve['macro_MAE_skill_vs_line_baseline_supported_cells']))
    assert max(errors)<1e-6,(max(errors),cells)
    result={'passed':True,'checkpoint_sha256':receipt['checkpoint_sha256'],'prediction_sha256':receipt['prediction_sha256'],
       'training_and_prediction_identities_exact':True,'training_prediction_groups_disjoint':True,'prediction_rows':len(te),
       'per_cell':cells,'independent_supported_cell_macro_raw_MAE_skill':macro,'maximum_metric_replay_error':max(errors),
       'raw_observed_magnitude_not_detection_or_latent_target':True,'measurement_SE_used_to_reweight_test_error':False,
       'new_study_outcomes_read':False,'independent_study_success_established':False,'candidate_promoted':False,
       'scope':'Saved prediction identity, hashes and all actual-magnitude MAE calculations. Neural forward replay is a separate checkpoint audit.',
       'code_sha256':sha(Path(__file__))}
    target=OUT/f'{stem}_epoch{args.epoch}_independent_saved_metric_audit.json';target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'macro_effect_size_skill':macro,'maximum_metric_replay_error':max(errors),'prediction_rows':len(te)}),flush=True)


if __name__=='__main__':main()
