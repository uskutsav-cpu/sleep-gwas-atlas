"""Freeze final settings once; sequence fits may run while public AVI inputs arrive."""
import sys,subprocess,time,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
DEST=ROOT/'outputs/sleep_frailty_ml_next';FINAL=DEST/'validation'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--arm',choices=['sequence','AVI'],required=True);args=parser.parse_args();FINAL.mkdir(exist_ok=True)
    assessment=json.loads((DEST/'coverage_development_assessment.json').read_text());assert assessment['all_fits_complete']
    plan=json.loads((DEST/'coverage_development_plan.json').read_text());work=Path(__file__).parent
    lock={'architecture':'Unchanged existing three-reader small response head; AVI adds only three zero-initialized hidden input columns. Author readers stay frozen.',
          'epochs':assessment['selected']['common_epoch'],'seeds':[42,7,2026],'held_studies':plan['held_studies'],
          'source_data_sha256':sha(DEST/'coverage_training_arrays.npz'),'pairs_sha256':sha(DEST/'coverage_measured_pairs.npz'),
          'selected_development_sha256':sha(DEST/'coverage_development_assessment.json'),'fit_code_sha256':{n:sha(work/n) for n in ['refit_coverage_candidate.py','matched_coverage_head.py','cell_matched_effect_network.py','coverage_inputs.py']},
          'no_further_source_architecture_or_epoch_search_before_test':True,'AVI_comparison':'Same rows, common shared initial parameters, actual batches, loss, update count and seeds. Inputs are exact-allele permissive composite AVI; not experimental truth.',
          'confidence_check':'For outer-held-source development checks, train error models only on other sources original calibration-region predictions from that same outer-source-omitted base model. Divide those whole regions into error fit and error scale. Do not use outer source to fit either stage.',
          'new_test_protocol':{'primary_tests':['Tardaguila2026 K562 independent-publication','Siraj2026 OL33 held-library macro across three known lines'],
             'model_required_relative_MAE_gain_vs_cell_and_line_medians':.05,'paired_95CI_gain_lower_above':0.,'magnitude_Spearman_minimum':.1,
             'AVI_required_relative_MAE_gain_vs_sequence':.01,'AVI_gain_95CI_lower_above':0.,'AVI_requires_both_primary_tests_to_pass':True,
             'nearby_ordering_and_difference_Spearman_max_drop':.01,'nearby_minimum_pairs':30,'nearby_minimum_1Mb_windows':20,
             'bootstrap_replicates':2000,'bootstrap_seed':20261108,'bootstrap_unit':'Fixed 1Mb genomic windows; not ancestry-specific LD blocks.',
             'nominal_interval_coverage':.9,'minimum_accepted_examples':30,'minimum_accepted_1Mb_windows':20,
             'unknown_blood_lines':'Report raw diagnostics separately; no authorized predictions or validated cell-specific differences.',
             'same_paper_OL27_reserved_cases':'Descriptive only, never substitute for independent publication.'},
          'reserved_outcomes_used_for_selection':False,'no_project_biological_rank_release_from_reporter_assay_only':True}
    path=FINAL/'final_settings_lock.json'
    if path.exists():assert json.loads(path.read_text())==lock
    else:write_json(path,lock)
    out=FINAL/'fits';out.mkdir(exist_ok=True);start=time.monotonic();names=[]
    for source in lock['held_studies']:
        for seed in lock['seeds']:
            name=f'fixed_{source}_{args.arm}_{seed}';names.append(name)
            if (out/(name+'_record.json')).exists():continue
            print('START FIXED COVERAGE',name,flush=True);command=[sys.executable,'-u',str(work/'refit_coverage_candidate.py'),'--held-study',source,'--seed',str(seed)]
            if args.arm=='AVI':command.append('--avi')
            with (out/(name+'.log')).open('w') as handle:subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT,check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2'})
            print('DONE FIXED COVERAGE',name,flush=True)
    write_json(FINAL/(args.arm+'_fixed_execution_receipt.json'),{'fits':names,'count':len(names),'settings_lock_sha256':sha(path),'elapsed_seconds':time.monotonic()-start,'reserved_performance_calculated':False})
    print('FIXED ARM COMPLETE',args.arm,flush=True)

if __name__=='__main__':main()
