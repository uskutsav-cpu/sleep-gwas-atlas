"""Finite sequential development comparisons; fresh effects stay closed."""
import subprocess, time
from v8_common import *


def main():
    plan = {
        'recipes': ['frozen', 'pooled_tune'],
        'source_omission_tests': ['none', 'Kosicki2025', 'IGVF240K', 'Salomon', 'Siraj'],
        'fixed_seed': 42, 'fixed_epochs': 6,
        'unseen_exact_cell_groups': 'Report raw diagnostics but exclude unsupported cells from the macro same-cell selection score; unfamiliar cells are not validated predictions.',
        'selected_recipe_rule': 'Average primary magnitude-error skill across whole-study omission tests with same-cell support, separately reporting random-locus validation and all local/cross-cell outcomes. Do not select from fresh Matoba results.',
        'additional_capacity_check': 'Two fixed end-to-end epochs on original training split, updating the last convolutional stage. This is a development experiment, not final testing.',
        'AlphaGenome_and_confidence_checks_remain_required': True,
        'prior_within_locus_pooled_tune_fit_already_observed': True,
        'documented_before_whole_study_omission_performance': True,
    }
    plan_path = OUT/'finite_development_comparison_plan.json'
    if not plan_path.exists():
        write_json(plan_path, plan)
    else:
        assert json.loads(plan_path.read_text()) == plan
    started = time.monotonic()
    for recipe in plan['recipes']:
        for study in plan['source_omission_tests']:
            name = f'{recipe}_{study}_42'
            record = OUT/(name+'_record.json')
            if record.exists():
                r = json.loads(record.read_text())
                assert r['completed_epochs'] == 6 and sha(r['checkpoint']) == r['checkpoint_sha256']
                print('EXISTING COMPLETED FIT VERIFIED', name, flush=True)
                continue
            log = OUT/(name+'.log')
            print('STARTING DEVELOPMENT COMPARISON', name, flush=True)
            with log.open('w') as f:
                subprocess.run([sys.executable, str(WORK/'train_development.py'), '--recipe', recipe, '--held-study', study, '--epochs', '6', '--device', 'cpu'], env=environment(), stdout=f, stderr=subprocess.STDOUT, check=True)
            r = json.loads(record.read_text())
            print('COMPLETED DEVELOPMENT COMPARISON', name, r['best_development_score'], 'elapsed', round(time.monotonic()-started, 1), flush=True)
    print('ALL FINITE POOLED DEVELOPMENT COMPARISONS COMPLETE', flush=True)
    # No study omission is claimed here: this checks a larger trainable
    # representation on untouched development loci, before fresh-study use.
    name = 'end_to_end_none_42'
    if not (OUT/(name+'_record.json')).exists():
        with (OUT/(name+'.log')).open('w') as f:
            subprocess.run([sys.executable, str(WORK/'train_development.py'), '--recipe', 'end_to_end', '--held-study', 'none', '--epochs', '2', '--device', 'mps'], env=environment(), stdout=f, stderr=subprocess.STDOUT, check=True)
    print('INITIAL V8 DEVELOPMENT CAPACITY CHECKS COMPLETE; FRESH VALIDATION STILL REQUIRED', flush=True)


if __name__ == '__main__':
    main()
