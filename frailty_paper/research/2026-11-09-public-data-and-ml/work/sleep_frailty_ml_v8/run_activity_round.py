"""Execute the declared finite sequence of measured-data comparisons."""
import subprocess, time
from v8_common import *


def main():
    plan = json.loads((OUT/'measured_activity_round_plan.json').read_text())
    start = time.monotonic()
    for recipe in plan['recipes']:
        for study in plan['held_out_development_studies']:
            name = f'round2_{recipe}_{study}_development_sequence_42'
            path = OUT/(name+'_record.json')
            if path.exists():
                r = json.loads(path.read_text())
                assert r['completed_epochs'] == plan['fixed_epochs'] and sha(r['checkpoint']) == r['checkpoint_sha256']
                print('VERIFIED COMPLETED ACTIVITY FIT', name, flush=True)
                continue
            print('STARTING DECLARED ACTIVITY FIT', name, flush=True)
            with (OUT/(name+'.log')).open('w') as f:
                subprocess.run([sys.executable, str(WORK/'train_activity_round.py'), '--recipe', recipe, '--held-study', study, '--epochs', str(plan['fixed_epochs'])], env=environment(), stdout=f, stderr=subprocess.STDOUT, check=True)
            r = json.loads(path.read_text()); assert sha(r['checkpoint']) == r['checkpoint_sha256']
            print('COMPLETED ACTIVITY FIT', name, 'total_seconds', round(time.monotonic()-start, 1), flush=True)
    write_json(OUT/'activity_round_execution_receipt.json', {'declared_plan_sha256': sha(OUT/'measured_activity_round_plan.json'), 'completed_fits': len(plan['recipes'])*len(plan['held_out_development_studies']), 'fresh_study_performance_read': False, 'elapsed_seconds': time.monotonic()-start})
    print('DECLARED ACTIVITY COMPARISONS COMPLETE; FRESH VALIDATION STILL CLOSED', flush=True)


if __name__ == '__main__':
    main()
