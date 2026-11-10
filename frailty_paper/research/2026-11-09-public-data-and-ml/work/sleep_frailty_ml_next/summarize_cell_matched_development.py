"""Record the finite comparison without turning development gains into validation."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *

DEST=ROOT/'outputs/sleep_frailty_ml_next'
SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'


def main():
    assert_previous_sources();execution=json.loads((DEST/'cell_matched_development_execution_receipt.json').read_text());assert execution['count']==16 and execution['matched_head_rows_and_initial_weights_verified']
    assert json.loads((DEST/'fitted_head_independent_audit.json').read_text())['passed'];plan=json.loads((DEST/'finite_cell_matched_development_plan.json').read_text());sources=plan['held_studies'][1:];comparisons=[];best_by_recipe={}
    for recipe in plan['recipes']:
        curves={s:json.loads((DEST/'development'/f'{recipe}_{s}_42_curve.json').read_text()) for s in sources}
        choices=[]
        for i in range(plan['epochs']):
            per={s:c[i]['macro_MAE_skill_vs_line_baseline_supported_cells'] for s,c in curves.items()};r={'recipe':recipe,'common_epoch':i+1,'mean_source_line_median_error_skill':float(np.mean(list(per.values()))),'per_source':per};choices.append(r);comparisons.append(r)
        best_by_recipe[recipe]=max(choices,key=lambda x:x['mean_source_line_median_error_skill'])
    selected=max(comparisons,key=lambda x:x['mean_source_line_median_error_skill']);pairs={}
    for recipe,r in best_by_recipe.items():pairs[recipe]=json.loads((DEST/'development'/f'{recipe}_none_42_curve.json').read_text())[r['common_epoch']-1]['comparisons']
    source=json.loads((SCREEN/'new_source_physical_and_historical_overlap_audit.json').read_text());pretrain=json.loads((SCREEN/'full_published_pretraining_exposure_audit.json').read_text());coverage=json.loads((SCREEN/'candidate_pretrained_readers/actual_pretraining_library_coverage_audit.json').read_text());libraries=json.loads((SCREEN/'Siraj_unused_libraries/shared_core_and_historical_overlap_audit.json').read_text())
    result={'selected_development_followup':selected,'best_by_recipe':best_by_recipe,'all_common_epoch_comparisons':comparisons,'held_locus_pair_diagnostics':pairs,
            'interpretation':'Only small, selected development gains. More epochs did not consistently help. Variant-specific cross-cell learning remains weak. New independent-study or held-library performance has not been calculated.',
            'new_sources':{'blood_after_physical_and_history_checks':source['accepted_tile_comparisons'],'blood_after_conservative_whole_published_pretraining_neighborhood_exclusion':pretrain['strict_no_exposure_designs'],
                           'blood_primary_candidate_distinct_variants':pretrain['strict_no_exposure_distinct_variants'],'actual_pretraining_DNA_covered_by_published_library_audit':coverage['all_actual_sequences_covered'],
                           'three_cell_library_roles_after_supervised_history_and_shared_DNA':libraries['by_library'],'new_library_eligibility_still_requires_pretraining_and_assay_audits':True},
            'pretrained_reader_weights_altered_by_this_training':False,'small_response_head_trained_by_backpropagation':True,
            'scientific_objective_complete':False,'Nature_Communications_genetic_manuscript_ready':False,'new_source_model_performance_used':False,
            'V9_adverse_external_results_or_AVI_rejection_changed':False,'next_steps':['Finish unused-library pretraining and coefficient/QC audits','At most one coverage expansion using eligible OL27 development regions, keeping OL33 and new blood tests outside training','Freeze model, uncertainty, matched AVI arms and predictions; then evaluate the new tests without another architecture/source search first']}
    write_json(DEST/'finite_development_assessment.json',result)
    table=[]
    for r in comparisons:
        table.append({'recipe':r['recipe'],'common_epoch':r['common_epoch'],'mean_source_line_median_error_skill':r['mean_source_line_median_error_skill'],**r['per_source']})
    write_csv(DEST/'whole_source_development_comparisons.csv',table)
    full=best_by_recipe['three_readers'];wtc=best_by_recipe['WTC11_only'];pair=pairs['three_readers']
    text=f'''CELL-MATCHED MODEL DEVELOPMENT: CURRENT STATUS

What changed
Three real author-trained DNA readers (blood, liver and stem cells) now process the full actual source DNA, preserving its length and assay orientation. No DNA was padded or extended. A small 8,461-parameter response head was trained by backpropagation. The pretrained readers themselves stayed fixed. The matched control masked the blood and liver features; training rows, initial head weights and update schedules were identical.

What was tested
Sixteen development fits: two input versions, eight epochs and seven whole-source omissions plus held loci. The source dataset remains the earlier revealed V9 measurements. No newly reserved blood or unused-library effect has been used to train or select these models. Computation and fitted gradients were checked independently.

Actual results
The best three-reader setup improved development error over the cell/culture median by {100*full['mean_source_line_median_error_skill']:.3f}% on average across seven omitted sources, at one common epoch ({full['common_epoch']}). The best stem-cell-only control improved it by {100*wtc['mean_source_line_median_error_skill']:.3f}% at its one common epoch ({wtc['common_epoch']}). These are selected development figures, not new independent validation. Further epochs did not consistently improve transfer.
On held development loci, nearby ordering accuracy was {100*pair['local']['ordering_accuracy_nonzero_measured_differences']:.2f}%. Within a fixed study/cell contrast, the average between-cell correlation was only {pair['cross']['macro_within_fixed_study_cell_pair_Spearman']:.3f}. Useful variant-specific cell differences remain unproven.

New free data verified
The independent-publication blood source has 354 verified paired designs after historical overlap checks, measured in four cell models. A conservative screen against the complete published pretraining libraries leaves 236 designs (49 distinct variants in 46 regions), potentially 944 measurements. All actual author pretraining sequences were found in the published libraries; unmapped/noncanonical genomic annotations still limit neighborhood claims. Original allele comparators and source filtered-count tables were checked; raw barcode processing and original GLM fitting have not been independently rebuilt.
Unused Siraj2026 libraries provide 13,817 OL27 candidate development designs and 18,219 OL33 candidate held-library designs after physical-reference, spatial-history and exact shared-DNA guards. Every design is measured in HepG2, K562 and SK-N-SH. Pretraining exposure and assay/target-availability checks remain. These are libraries from the same paper and lab, not independent publications or donors.

Bounded follow-up
Finish those source checks, make at most one declared coverage expansion with the eligible OL27 development regions, then freeze the model, confidence rules, both matched AVI arms and all predictions before evaluating the new blood study and held OL33 library. Do not make another architecture/source change before that evaluation. A weak result must be retained as a weak result.

Still unresolved
The prior V9 model failed its untouched accuracy and confidence tests. Its frozen results and AVI rejection remain unchanged. The new models have not passed an external test or solved confidence across experiments. No project biological ranking is authorized. This is laboratory molecular research, not a predictor of insomnia/frailty in people, independent genetic-score validation, causal-gene evidence or a publication-readiness claim.
'''
    (DEST/'PLAIN_LANGUAGE_PROGRESS.txt').write_text(text);write_json(DEST/'bounded_followup_plan.json',{'maximum_additional_coverage_expansions_before_new_external_test':1,'source':'Eligible OL27 development components only; all OL33-connected and blood-test components excluded.',
          'architecture_search_before_next_test':'No additional search after the coverage-expanded candidate is fixed. Retain one existing response-head architecture and a finite common-epoch selection.',
          'new_test_outcomes_used':False,'required_before_test':['Finish all source/pretraining/QC checks','Freeze model, uncertainty, matched AVI settings and both arms predictions'],'scope':'Do not count same-paper libraries as independent publications; do not interpret zero answers as useful prediction.'})
    print(text,flush=True)


if __name__=='__main__':main()
