"""Summarize all fixed development folds without evaluating reserved data.

Separate overall assay detection from ranking individual variants within active
regions. Do not pool predictions from differently fitted folds for selection.
"""
import argparse,json,csv,subprocess,numpy as np
from common import *
FOLDS=[list(range(4,10)),list(range(10,16)),list(range(16,23))]
def families(val,rep):
    suffix='-'.join(map(str,val));return [
        ('contrast_full_val'+suffix,{'_ensemble_effect':'contrast_no_avi','_ensemble_effect_avi':'contrast_avi_teacher'}),
        ('motif_chromatin_full_val'+suffix,{'_ensemble_activity_effect':'motif_chromatin_no_avi','_ensemble_activity_avi':'motif_chromatin_avi_teacher'}),
        (rep+'_full_avi_input_val'+suffix,{'_ensemble_activity_effect':'biological_avi_input'}),
        ('development_baselines_val'+suffix,None)]
def main():
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');a=p.parse_args()
    if not a.partial:assert (OUT/'development_robustness_finished.json').is_file(),'Fitting must finish before final comparison'
    plan=json.loads((OUT/'development_robustness_locked_config.json').read_text());rep=plan['screening_selected_AVI_input_representation'];data,_,manifest,_,_,_=load();active=np.array([int(r['active_screen']) for r in manifest]);fold_results=[];rchecks=[]
    for val in FOLDS:
        for prefix,mapping in families(val,rep):
            if not (OUT/(prefix+'_selection.json')).is_file():continue
            rows=json.loads((OUT/(prefix+'_search_results.json')).read_text());lookup={r['name']:r for r in rows}
            with np.load(OUT/(prefix+'_validation_predictions.npz')) as z:
                ix=z['row_index'];assert set(data['chromosome'][ix])==set(val) and np.all(data['chromosome'][ix]>=4)
                assert np.array_equal(z['label'],data['label'][ix]);preds={k:z[k] for k in lookup}
            conditional={}
            cp=OUT/(prefix+'_conditional_validation_predictions.npz')
            if cp.exists():
                with np.load(cp) as z:assert np.array_equal(z['row_index'],ix);conditional={k:z[k] for k in lookup}
            names={name:(prefix+name if prefix+name in lookup else None) for name in (mapping or {})}
            if mapping:chosen=[(mapping[suffix],name) for suffix,name in names.items() if name]
            else:chosen=[('baseline_'+r['screening_source_record'],r['name']) for r in rows]
            for family,name in chosen:
                overall=preds[name];variant=conditional.get(name,overall);allmetrics=study_metrics(data['label'][ix],overall,data['study'][ix]);within={}
                for s,study in enumerate(STUDIES):
                    mask=(data['study'][ix]==s)&(active[ix]==1);within[study]=metrics(data['label'][ix[mask]],variant[mask])
                    # Independent R replay includes full and active-only predictions.
                    for stratum,use,score in [('all',data['study'][ix]==s,overall),('active',mask,variant)]:
                        path=OUT/(f'development_R_{"-".join(map(str,val))}_{family}_{s}_{stratum}.csv')
                        with path.open('w',newline='') as f:
                            w=csv.writer(f);w.writerow(['row_id','label','score']);w.writerows(zip(ix[use],data['label'][ix[use]].astype('int32'),score[use]))
                        target=path.with_name(path.stem+'_metrics.csv');subprocess.run(['Rscript',str(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R'),str(path),str(target)],capture_output=True,text=True,check=True);r=list(csv.DictReader(target.open()))[0];expected=allmetrics['studies'][study] if stratum=='all' else within[study];errors={k:abs(float(r[k])-expected[k]) for k in ['auroc','average_precision']};assert max(errors.values())<1e-10,(family,study,stratum,errors);rchecks.append({'fold':val,'family':family,'study':study,'stratum':stratum,'errors':errors})
                fold_results.append({'fold':val,'family':family,'record_name':name,'all_calls':allmetrics,'conditional_active':within,'macro_conditional_AP_skill':float(np.mean([m['AP_skill_over_prevalence'] for m in within.values()])),'variant_score_definition':'Active-region conditional expert' if name in conditional else 'Overall call score restricted to measured-active evaluation rows','members':lookup[name]['config'].get('members',[]),'artifact':lookup[name]['artifact']})
    aggregate=[]
    for family in sorted({r['family'] for r in fold_results}):
        rows=[r for r in fold_results if r['family']==family]
        aggregate.append({'family':family,'completed_folds':len(rows),'all_fold_mean_all_call_AP_skill':float(np.mean([r['all_calls']['macro_AP_skill'] for r in rows])),'all_fold_mean_conditional_AP_skill':float(np.mean([r['macro_conditional_AP_skill'] for r in rows])),'per_study_mean':{s:{'all_call_auroc':float(np.mean([r['all_calls']['studies'][s]['auroc'] for r in rows])),'all_call_AP':float(np.mean([r['all_calls']['studies'][s]['average_precision'] for r in rows])),'conditional_auroc':float(np.mean([r['conditional_active'][s]['auroc'] for r in rows])),'conditional_AP':float(np.mean([r['conditional_active'][s]['average_precision'] for r in rows]))} for s in STUDIES}})
    complete=all(r['completed_folds']==3 for r in aggregate) and len(aggregate)==9
    if not a.partial:assert complete,(len(aggregate),[(r['family'],r['completed_folds']) for r in aggregate])
    output={'status':'Complete development-only comparison' if complete else 'Incomplete development preview; no final family selection','fold_results':fold_results,'aggregate':aggregate,'R_replays':rchecks,'test_evaluated':False,'qualification':'Folds served for model development, early stopping and family selection after4–5 screening. These are not unbiased nested-CV estimates or independent biological replication. Reserved Siraj chromosomes1–3 remain necessary. Each study is reported separately.'}
    path=OUT/('development_comparison_partial.json' if a.partial else 'development_comparison.json');path.write_text(json.dumps(output,indent=2)+'\n');print('Development comparisons saved',len(fold_results),'fold/family results;',len(rchecks),'independent R metric replays',flush=True)
if __name__=='__main__':main()
