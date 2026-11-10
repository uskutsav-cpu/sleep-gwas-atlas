"""Execute the original finite comparison through its attested read-only fitter.

The original role list, matching checks, epoch selection and stopping gates are
imported unchanged. This filename preserves the shared driver's live-PID lock.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PARENT))
import run_resumable_comparison as original
from process_control import active_fit, command, driver_lock
from readonly_native_runtime import execution_revision
from resumable_common import DEST, sha


def fit(recipe,source,plan):
    directory = DEST/'development'
    directory.mkdir(exist_ok=True)
    name = f'{recipe}_{source}_42'
    if (directory/(name+'_record.json')).exists():
        return
    script = str(PARENT/'readonly_native_runtime.py')
    old_script = str(PARENT/'train_resumable_development.py')
    workers = [(path,active_fit(path,recipe,source)) for path in [script,old_script]]
    workers = [(path,pid) for path,pid in workers if pid]
    assert len(workers) <= 1, workers
    if workers:
        path,pid = workers[0]
        print('FOLLOW VERIFIED EXISTING FIT',name,pid,flush=True)
        while path in command(pid):
            time.sleep(10)
        if (directory/(name+'_record.json')).exists():
            return
    started = (directory/(name+'_started.json')).exists()
    if started:
        assert (directory/(name+'_recovery.pt')).exists(), name
    print('START ATTESTED READ-ONLY NATIVE FIT',name,'resume',started,flush=True)
    with (directory/(name+'.log')).open('a' if started else 'w') as handle:
        subprocess.run([sys.executable,'-u',script,'--recipe',recipe,'--held-study',source]+(['--resume'] if started else []),
                       check=True,stdout=handle,stderr=subprocess.STDOUT,
                       env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2',
                            'OPENBLAS_NUM_THREADS':'2','VECLIB_MAXIMUM_THREADS':'2',
                            'PYTORCH_MPS_HIGH_WATERMARK_RATIO':'0.4','PYTORCH_MPS_LOW_WATERMARK_RATIO':'0.3'})
    print('DONE ATTESTED READ-ONLY NATIVE FIT',name,flush=True)


def main():
    _,revision_hash = execution_revision()
    proof = json.loads((DEST/'readonly_runtime_recovery_audit.json').read_text())
    assert proof['passed'] and proof['execution_revision_sha256'] == revision_hash
    original.fit = fit
    receipt = {'PID':os.getpid(),'driver_source_sha256':sha(Path(__file__)),
               'inherited_finite_driver_source_sha256':sha(Path(original.__file__)),
               'execution_revision_sha256':revision_hash,
               'original_roles_matching_epoch_selection_and_gates_imported_unchanged':True,
               'comparison_restarted':False,'scientific_goal_complete':False}
    (DEST/'readonly_finite_driver_execution.json').write_text(json.dumps(receipt,indent=2)+'\n')
    original.main()


if __name__ == '__main__':
    with driver_lock(DEST/'driver_lock.json'):
        main()
