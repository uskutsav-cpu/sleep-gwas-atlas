"""Select on development only, then refit fixed recipes without final outcomes."""
import sys,argparse,subprocess
from train_development import *
def select(d,pairs):
    assert not (OUT/'final_test_evaluation_started.json').exists()
    records=json.loads((OUT/'development_records.json').read_text());assert len(records)==17
    va=np.flatnonzero(d['split']==1);families={};predictions={}
    for name in sorted({r['name'] for r in records}):
        members=[r for r in records if r['name']==name];parts=[]
        for r in members:
            assert sha(r['checkpoint'])==r['checkpoint_sha256'];path=Path(r['checkpoint']);path=path.with_name(path.stem+'_validation.npz')
            with np.load(path) as z:assert np.array_equal(z['row_index'],va);parts.append(z['prediction'])
        p=np.mean(parts,axis=0);predictions[name]=p;families[name]={'members':members,'metrics':evaluate(d,va,p,pairs,name=='contrast_neural_only'),'primary_eligible':name!='contrast_neural_only'}
    eligible=sorted([name for name in families if families[name]['primary_eligible']],key=lambda n:families[n]['metrics']['selection_score'],reverse=True);candidates=[]
    for name in eligible:candidates.append({'name':name,'family_weights':{name:1.},'metrics':families[name]['metrics']})
    for weight in [.5,.75]:
        a,b=eligible[:2];p=weight*predictions[a]+(1-weight)*predictions[b];name=f'mixture_top2_{int(weight*100)}';predictions[name]=p;candidates.append({'name':name,'family_weights':{a:weight,b:1-weight},'metrics':evaluate(d,va,p,pairs)})
    candidates.sort(key=lambda r:r['metrics']['selection_score'],reverse=True);winner=candidates[0]
    report={'selected':winner,'candidates':candidates,'families':families,'validation_row_indices_sha256':hashlib.sha256(va.tobytes()).hexdigest(),'selection_plan_sha256':sha(OUT/'model_search_plan.json'),'final_outcomes_opened':False,'scope':'Development selection. No final-study result used.'};write_json(OUT/'development_selection.json',report)
    np.savez_compressed(OUT/'development_ensemble_predictions.npz',row_index=va,**predictions)
    return report
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-for-training-pid',type=int);parser.add_argument('--select-only',action='store_true');args=parser.parse_args()
    if args.wait_for_training_pid:
        while True:
            r=subprocess.run(['ps','-p',str(args.wait_for_training_pid),'-o','command='],capture_output=True,text=True)
            if r.returncode or 'sleep_frailty_ml_v6/train_development.py' not in r.stdout:break
            time.sleep(10)
    d,receipt,bank_x,pairs=prepare_runtime();report=select(d,pairs);print('V5 SELECTED DEVELOPMENT RECIPE',report['selected']['name'],report['selected']['metrics']['selection_score'],flush=True)
    if args.select_only:return
    mx.metal.set_memory_limit(1536*1024**2);mx.metal.set_cache_limit(64*1024**2);plan=freeze_plan();primary=set(report['selected']['family_weights']);controls={'local_mlp','contrast_standard','contrast_hard','contrast_neural_only'};needed=primary|controls;records=[]
    for name in sorted(needed):
        family=report['families'][name]
        for old in family['members']:
            if name.startswith('tree_'):path=OUT/f'final_{name}_record.json'
            else:path=OUT/f'final_{name}_{old["seed"]}_record.json'
            if path.exists():r=json.loads(path.read_text());assert sha(r['checkpoint'])==r['checkpoint_sha256']
            elif name.startswith('tree_'):r=fit_tree(d,receipt,bank_x,pairs,old['depth'],stage='final')
            else:r=fit_neural(d,receipt,bank_x,pairs,old['recipe'],old['seed'],plan,stage='final',fixed_epochs=old['selected_epoch'])
            records.append(r);write_json(OUT/'final_fit_records.json',records);mx.metal.clear_cache()
    registry={'selected_family_weights':report['selected']['family_weights'],'selection_sha256':sha(OUT/'development_selection.json'),'feature_receipt_sha256':sha(OUT/'fixed_features_receipt.json'),'records':records,'calibration_rows_fitted':False,'final_study_fitted_or_evaluated':False,'primary_AlphaGenome_inputs_or_weights':False,'scope':'Refitted DNA-effect laboratory models. Separate probability/uncertainty calibration and final-study evaluation still required.'}
    write_json(OUT/'final_fit_registry.json',registry);print('FINAL REFITS COMPLETE',len(records),flush=True)
if __name__=='__main__':main()
