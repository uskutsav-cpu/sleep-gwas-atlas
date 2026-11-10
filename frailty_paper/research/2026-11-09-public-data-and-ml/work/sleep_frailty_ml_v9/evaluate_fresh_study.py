"""One frozen-model evaluation; no threshold or model selection on fresh effects."""
from final_validation_common import *
from scipy.stats import spearmanr
import importlib.util


def rho(y,p):
    return float(spearmanr(y,p).statistic) if len(y)>2 and np.ptp(y)>0 and np.ptp(p)>0 else None


def interval(values):
    v=np.asarray([x for x in values if x is not None and np.isfinite(x)])
    return np.quantile(v,[.025,.975]).tolist() if len(v) else [None,None]


def paired_bootstrap(y,p,baseline,cell_baseline,groups,mask=None):
    if mask is None:mask=np.ones(len(y),dtype=bool)
    y,p,baseline,cell_baseline,groups=[a[mask] for a in [y,p,baseline,cell_baseline,groups]]
    if not len(y):return {'n':0,'windows':0}
    error=abs(y-p);be=abs(y-baseline);ce=abs(y-cell_baseline);unique=np.unique(groups);blocks=[np.flatnonzero(groups==s) for s in unique]
    rng=np.random.default_rng(20261108);line_gain=[];cell_gain=[];rhos=[];relative=[]
    for _ in range(2000):
        ix=np.concatenate([blocks[j] for j in rng.integers(len(blocks),size=len(blocks))]);line_gain.append(float((be[ix]-error[ix]).mean()));cell_gain.append(float((ce[ix]-error[ix]).mean()));relative.append(float(1-error[ix].mean()/be[ix].mean()));rhos.append(rho(y[ix],p[ix]))
    return {'n':len(y),'windows':len(unique),'MAE':float(error.mean()),'line_or_cell_baseline_MAE':float(be.mean()),'cell_baseline_MAE':float(ce.mean()),
            'relative_MAE_improvement_vs_line_baseline':float(1-error.mean()/be.mean()),'relative_MAE_improvement_vs_cell_baseline':float(1-error.mean()/ce.mean()),
            'absolute_error_improvement_vs_line_baseline_CI95':interval(line_gain),'absolute_error_improvement_vs_cell_baseline_CI95':interval(cell_gain),
            'relative_MAE_improvement_vs_line_baseline_CI95':interval(relative),'magnitude_Spearman':rho(y,p),'magnitude_Spearman_CI95':interval(rhos),
            'bootstrap':'2,000 paired 1-Mb-window resamples; windows approximate dependence and are not ancestry-specific LD blocks.'}


def nearby_pairs(rows):
    coords=[(r['variant_GRCh38'].split(':')[0],int(r['variant_GRCh38'].split(':')[1])) for r in rows];order=sorted(range(len(rows)),key=lambda j:coords[j]);pairs=[]
    for a,b in zip(order[:-1],order[1:]):
        if coords[a][0]==coords[b][0] and 0<coords[b][1]-coords[a][1]<=25:pairs.append([a,b])
    return np.asarray(pairs,dtype='int64').reshape(-1,2)


def local_metrics(y,p,pairs):
    if not len(pairs):return {'pairs':0,'ordering_accuracy':None,'difference_Spearman':None}
    a,b=pairs.T;dy=y[a]-y[b];dp=p[a]-p[b];ok=dy!=0
    return {'pairs':len(pairs),'nonzero_measured_difference_pairs':int(ok.sum()),'difference_MAE':float(abs(dy-dp).mean()),'zero_difference_MAE':float(abs(dy).mean()),
            'ordering_accuracy':float(np.mean(np.sign(dy[ok])==np.sign(dp[ok]))) if ok.any() else None,'difference_Spearman':rho(dy,dp)}


def primary_gate(result):
    return bool(result['relative_MAE_improvement_vs_line_baseline']>=.05 and result['relative_MAE_improvement_vs_cell_baseline']>=.05
                and result['absolute_error_improvement_vs_line_baseline_CI95'][0]>0 and result['absolute_error_improvement_vs_cell_baseline_CI95'][0]>0
                and result['magnitude_Spearman'] is not None and result['magnitude_Spearman']>=.1 and result['magnitude_Spearman_CI95'][0]>0)


def main():
    assert_previous_sources();lockpath=FINAL/'before_fresh_outcomes_lock.json';assert lockpath.exists();assert not (FINAL/'fresh_outcome_evaluation.json').exists()
    lock=json.loads(lockpath.read_text())
    for name,digest in lock['files'].items():assert sha(name)==digest,name
    protocol=json.loads((FINAL/'external_validation_protocol.json').read_text());source=V8/'sources/Matoba_journal_Supplementary_Table3.xlsx'
    assert sha(source)==protocol['official_journal_source_sha256']
    spec=importlib.util.spec_from_file_location('verified_readonly_source',ROOT/'work/sleep_frailty_ml_v7/extract_neuroimmune_sources.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    table=reader.read_tables(source,['b. Allelic effect (Veh)'])['b. Allelic effect (Veh)']['rows'];header=table[0]
    assert 'SNP' in header and 'logFC' in header,header
    effects={}
    for row in table[1:]:
        r=dict(zip(header,row));assert r['SNP'] not in effects;effects[r['SNP']]=abs(float(r['logFC']))
    allrows=list(csv_rows(V8/'blind_inference_row_metadata.csv'));fresh=[r for r in allrows if r['role']=='fresh_test'];assert len(fresh)==610
    y=np.array([effects[r['rsid']] for r in fresh]);assert np.isfinite(y).all()
    # Outcome extraction happens only after every model, interval and acceptance
    # decision has been committed by hash in the pre-outcome lock.
    write_json(FINAL/'fresh_outcome_join_receipt.json',{'lock_sha256':sha(lockpath),'journal_source_sha256':sha(source),'all_661_source_coefficients_considered':len(effects)==661,
               'eligible_coefficients':len(y),'target':'Absolute original fitted log2 allelic coefficient; signed allele convention not inferred.',
               'outcome_sha256':hashlib.sha256(y.tobytes()).hexdigest(),'outcome_or_significance_filtering':False})
    groups=np.array([r['variant_GRCh38'].split(':')[0]+':'+str((int(r['variant_GRCh38'].split(':')[1])-1)//1000000) for r in fresh]);pairs=nearby_pairs(fresh)
    local_windows=len(np.unique(groups[pairs[:,0]])) if len(pairs) else 0;arms={};predictions={}
    for arm in ['sequence','AVI']:
        with np.load(FINAL/f'{arm}_blinded_predictions.npz') as z:raw={k:z[k] for k in z.files}
        ids=np.array([j for j,i in enumerate(raw['row_index']) if allrows[i]['role']=='fresh_test']);assert [allrows[raw['row_index'][j]]['rsid'] for j in ids]==[r['rsid'] for r in fresh]
        pred={k:v[ids] for k,v in raw.items() if k!='row_index'};predictions[arm]=pred;p=pred['magnitude'];hw=pred['interval_halfwidth'];covered=abs(y-p)<=hw;accepted=pred['pretest_policy_accepted'].astype(bool)
        metrics=paired_bootstrap(y,p,pred['baseline'],pred['cell_baseline'],groups);selected=paired_bootstrap(y,p,pred['baseline'],pred['cell_baseline'],groups,accepted)
        qualified=primary_gate(metrics);coverage=float(covered.mean());selected_coverage=float(covered[accepted].mean()) if accepted.any() else None
        selected_positive=bool(selected.get('magnitude_Spearman') is not None and selected['magnitude_Spearman']>0)
        confidence_gate=bool(qualified and accepted.sum()>=30 and len(np.unique(groups[accepted]))>=20 and coverage>=.9 and selected_coverage is not None and selected_coverage>=.9
                             and selected.get('absolute_error_improvement_vs_line_baseline_CI95',[None])[0] is not None and selected['absolute_error_improvement_vs_line_baseline_CI95'][0]>0
                             and selected['absolute_error_improvement_vs_cell_baseline_CI95'][0]>0 and selected_positive)
        local=local_metrics(y,p,pairs)
        if len(pairs):
            a,b=pairs.T;paircal=json.loads((FINAL/f'{arm}_confidence_calibration.json').read_text())['paired']['local'];key='neural_progenitor|neural_progenitor';factor=paircal['by_cell_pair'].get(key,{'factor':paircal['global_factor']})['factor']
            halfwidth=factor*(pred['raw_error_scale'][a]+pred['raw_error_scale'][b]);dy=y[a]-y[b];dp=p[a]-p[b];answers=accepted[a]&accepted[b]&(abs(dp)>halfwidth)
            local.update(direct_pair_interval_coverage=float(np.mean(abs(dy-dp)<=halfwidth)),direct_pair_calibration_factor=factor,pretest_answerable_pairs=int(answers.sum()),authorized_pairs=int(answers.sum()) if confidence_gate else 0,
                         individual_intervals_used_to_claim_paired_guarantee=False)
        arms[arm]={'all_cases':metrics,'pretest_accepted_cases':selected,'all_case_interval_coverage':coverage,'accepted_case_interval_coverage':selected_coverage,
                   'mean_interval_width':float((p+hw-np.maximum(p-hw,0)).mean()),'pretest_policy_accepted_cases':int(accepted.sum()),'authorized_cases':int(accepted.sum()) if confidence_gate else 0,
                   'full_model_quality_gate_passed':qualified,'uncertainty_and_delivery_gate_passed':confidence_gate,'nearby_comparisons':local}
    seq,avi=predictions['sequence'],predictions['AVI'];seqerr=abs(y-seq['magnitude']);avierr=abs(y-avi['magnitude']);rng=np.random.default_rng(20261108);unique=np.unique(groups);blocks=[np.flatnonzero(groups==s) for s in unique];gains=[]
    for _ in range(2000):
        ix=np.concatenate([blocks[j] for j in rng.integers(len(blocks),size=len(blocks))]);gains.append(float((seqerr[ix]-avierr[ix]).mean()))
    improvement=float(1-avierr.mean()/seqerr.mean());ci=interval(gains);minimum=protocol['V9_minimum_local_confirmation'];local_enough=len(pairs)>=minimum['pairs'] and local_windows>=minimum['independent_1Mb_windows']
    sl,al=[arms[a]['nearby_comparisons'] for a in ['sequence','AVI']]
    no_order_drop=sl['ordering_accuracy'] is not None and al['ordering_accuracy'] is not None and al['ordering_accuracy']>=sl['ordering_accuracy']-.01
    no_rho_drop=sl['difference_Spearman'] is not None and al['difference_Spearman'] is not None and al['difference_Spearman']>=sl['difference_Spearman']-.01
    retain=bool(improvement>=.01 and ci[0]>0 and local_enough and no_order_drop and no_rho_drop)
    result={'version':'V9 finite, locked whole-study molecular validation','before_outcome_lock_sha256':sha(lockpath),'fresh_cases':len(y),'independent_1Mb_windows':len(unique),'arms':arms,
            'AVI_comparison':{'relative_MAE_improvement_vs_sequence':improvement,'paired_absolute_error_improvement_CI95':ci,'nearby_pairs':len(pairs),'nearby_1Mb_windows':local_windows,'local_confirmation_sufficient':bool(local_enough),
                              'ordering_drop_within_limit':bool(no_order_drop),'difference_rho_drop_within_limit':bool(no_rho_drop),'retain_AVI':retain,'attribution':'Composite AVI includes other predictors. No AlphaGenome-specific causal attribution.'},
            'scope':'Reporter-effect magnitude in the eligible assay-active primary-neural-progenitor subset; one donor, no established independence from the earlier primary-cell donor.',
            'not_established':['Independent human insomnia/frailty score validation','Disease causation','Correct signed original allele convention','Independent validation of differences across cell types','Applicability to unmeasured native project contexts','Full pretrained genomic-neighborhood independence'],
            'fresh_study_now_revealed':True,'future_model_change_requires_another_untouched_experiment':True}
    write_json(FINAL/'fresh_outcome_evaluation.json',result)
    output=[]
    for i,r in enumerate(fresh):
        row={**r,'measured_absolute_original_log2_effect':float(y[i]),'bootstrap_1Mb_window':groups[i]}
        for arm,p in predictions.items():
            for k in ['magnitude','baseline','cell_baseline','interval_halfwidth','pretest_policy_accepted']:row[f'{arm}_{k}']=float(p[k][i]) if k!='pretest_policy_accepted' else bool(p[k][i])
            row[f'{arm}_authorized_answer']=bool(p['pretest_policy_accepted'][i] and arms[arm]['uncertainty_and_delivery_gate_passed'])
        output.append(row)
    write_csv(FINAL/'fresh_study_predictions_and_observations.csv',output)
    assert_previous_sources();print(json.dumps(result),flush=True)


if __name__=='__main__':main()
