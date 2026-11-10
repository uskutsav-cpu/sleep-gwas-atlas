"""Matched refits and source-omitted ensembles; no fresh outcome evaluation."""
import subprocess, time
from v8_common import *


def main():
    lock = json.loads((OUT/'activity_round_final_settings_lock.json').read_text())
    completed = []; started = time.monotonic()
    for study in ['none']+lock['source_omission_studies']:
        for seed in lock['seeds']:
            matched = []
            for arm in lock['matched_input_arms']:
                name = f'fixed_{lock["selected_family"]}_{lock["selected_recipe"]}_{study}_{arm}_{seed}'
                path = OUT/(name+'_record.json')
                if path.exists():
                    record = json.loads(path.read_text())
                    assert sha(record['checkpoint']) == record['checkpoint_sha256'] and sha(record['predictions']) == record['predictions_sha256']
                    print('VERIFIED COMPLETED FIXED FIT', name, flush=True)
                else:
                    command = [sys.executable, str(WORK/'refit_fixed_candidate.py'), '--held-study', study, '--seed', str(seed)]
                    if arm == 'AVI': command.append('--avi')
                    print('STARTING MATCHED FIXED FIT', name, flush=True)
                    with (OUT/(name+'.log')).open('w') as f:
                        subprocess.run(command, env=environment(), stdout=f, stderr=subprocess.STDOUT, check=True)
                    record = json.loads(path.read_text())
                matched.append(record); completed.append(path.name)
            a, b = matched
            for key in ['shared_initial_parameter_sha256', 'training_indices_sha256', 'predict_indices_sha256', 'epoch_count', 'seed']:
                assert a[key] == b[key], (study, seed, key)
            print('MATCHED PAIR VERIFIED', study, seed, 'seconds', round(time.monotonic()-started, 1), flush=True)
    write_json(OUT/'fixed_refit_execution_receipt.json', {'fits': completed, 'completed_fits': len(completed), 'all_matched_initial_shared_weights_and_training_rows_verified': True, 'settings_lock_sha256': sha(OUT/'activity_round_final_settings_lock.json'), 'fresh_performance_calculated': False, 'elapsed_seconds': time.monotonic()-started})
    print('MATCHED FINAL AND SOURCE-OMITTED ENSEMBLES COMPLETE', flush=True)


if __name__ == '__main__':
    main()
