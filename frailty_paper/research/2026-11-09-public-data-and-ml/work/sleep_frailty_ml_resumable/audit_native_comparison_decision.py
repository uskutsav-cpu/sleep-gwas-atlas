"""Replay the locked decision using independently checked saved predictions."""
import datetime,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ml_resumable'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 plan_path=OUT/'finite_stable_native_plan.json';plan=read(plan_path)
 assessment_path=OUT/'native_comparison_assessment.json';assessment=read(assessment_path)
 evidence={str(plan_path.relative_to(ROOT)):sha(plan_path),str(assessment_path.relative_to(ROOT)):sha(assessment_path)}
 sources=plan['phase_one_roles'][1:];roles=plan['phase_one_roles'];recipes=plan['recipes']
 assert sources==['Salomon','Kosicki2025','Retallick2026'] and recipes==['frozen_reader','adapt_last_stage'] and plan['epochs']==3
 metrics={}
 for recipe in recipes:
  for source in roles:
   stem=f'{recipe}_{source}_42';record_path=OUT/'development'/(stem+'_record.json');record=read(record_path)
   assert record['completed_epochs']==3
   evidence[str(record_path.relative_to(ROOT))]=sha(record_path)
   metrics[recipe,source]=[]
   for epoch in range(1,4):
    path=OUT/f'{stem}_epoch{epoch}_independent_saved_metric_audit.json';row=read(path)
    assert row['passed'] and row['training_and_prediction_identities_exact'] and row['training_prediction_groups_disjoint']
    receipt=read(OUT/'development'/f'{stem}_epoch{epoch}_receipt.json')
    assert row['checkpoint_sha256']==receipt['checkpoint_sha256'] and row['prediction_sha256']==receipt['prediction_sha256']
    evidence[str(path.relative_to(ROOT))]=sha(path);metrics[recipe,source].append(row['independent_supported_cell_macro_raw_MAE_skill'])
 for source in roles:
  rows=[read(OUT/'development'/f'{recipe}_{source}_42_record.json') for recipe in recipes]
  for key in ['training_index_sha256','prediction_index_sha256','initial_parameter_sha256','normalization','context_transformations','actual_batch_schedule_sha256','actual_backprop_updates','completed_epochs','new_original_Siraj_signed_training_rows','new_original_neuronal_signed_training_rows']:
   assert rows[0][key]==rows[1][key],(source,key)
 choices=[]
 for epoch in range(1,4):
  a={s:metrics['adapt_last_stage',s][epoch-1] for s in sources};b={s:metrics['frozen_reader',s][epoch-1] for s in sources}
  choices.append({'epoch':epoch,'adapted_mean_skill':math.fsum(a.values())/3,'frozen_mean_skill':math.fsum(b.values())/3,'adapted_per_source':a,'frozen_per_source':b})
 selected=max(choices,key=lambda row:row['adapted_mean_skill']);epoch=selected['epoch'];pairs={}
 for recipe in recipes:
  path=OUT/f'{recipe}_none_42_epoch{epoch}_independent_pair_audit.json';row=read(path);assert row['passed'];pairs[recipe]=row['comparisons'];evidence[str(path.relative_to(ROOT))]=sha(path)
 a,b=pairs['adapt_last_stage'],pairs['frozen_reader']
 flags={'mean_actual_magnitude_skill_at_least_one_percent':selected['adapted_mean_skill']>=.01,'improvement_vs_matched_frozen_at_least_half_percentage_point':selected['adapted_mean_skill']-selected['frozen_mean_skill']>=.005,'no_neural_source_error_worsens_over_five_percent':all(1-selected['adapted_per_source'][s]<=1.05*(1-selected['frozen_per_source'][s]) for s in sources),'nearby_ordering_drop_at_most_point_zero_one':a['local']['ordering_accuracy_nonzero_measured_differences']>=b['local']['ordering_accuracy_nonzero_measured_differences']-.01,'within_fixed_cell_pair_rho_drop_at_most_point_zero_one':a['cross']['macro_within_fixed_study_cell_pair_Spearman']>=b['cross']['macro_within_fixed_study_cell_pair_Spearman']-.01,'nearby_difference_error_beats_zero_difference':a['local']['difference_MAE']<a['local']['zero_difference_MAE']}
 gate=assessment['gate'];assert epoch==gate['selected']['common_epoch'] and flags==gate['criteria'] and all(flags.values())==gate['passed']
 assert abs(selected['adapted_mean_skill']-gate['selected']['adapted_mean_neural_source_skill'])<1e-6
 assert abs(selected['frozen_mean_skill']-gate['selected']['frozen_mean_neural_source_skill'])<1e-6
 assert not gate['passed'] and not assessment['phase_two_run']
 result={'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_eight_fits_complete_and_matched':True,'all_24_saved_magnitude_results_independently_checked':True,'selected_epoch_pair_results_independently_checked':True,'selected_common_epoch':epoch,'independent_selected_results':selected,'independently_replayed_criteria':flags,'scientific_gate_passed':False,'original_stopping_decision_confirmed':True,'other_four_source_fits_not_authorized_by_this_comparison':True,'external_study_outcomes_used':False,'candidate_promoted':False,'goal_complete':False,'scope':'Independent saved-prediction metric and fixed-decision replay. Computational checks passed; scientific accuracy criteria failed. No clinical, causal, external-study or useful-confidence validation.','evidence_sha256':evidence,'code_sha256':sha(Path(__file__))}
 (OUT/'native_comparison_independent_decision_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['evidence_sha256','code_sha256']},indent=2))
if __name__=='__main__':main()
