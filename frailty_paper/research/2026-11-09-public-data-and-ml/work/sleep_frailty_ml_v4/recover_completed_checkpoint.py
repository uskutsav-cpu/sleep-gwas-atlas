"""Audit the completed checkpoint after the native run was interrupted.

This does not fit, resume, or select another model. The incomplete second
native epoch is explicitly excluded; only the previously saved checkpoint
and its recorded validation metrics are eligible for the assay probe.
"""
import subprocess,sys
from common import *

def main():
    live=subprocess.run(['ps','ax','-o','command='],capture_output=True,text=True,check=True).stdout.splitlines()
    assert not any(line.startswith(('python3 ', '/usr/local/bin/python3 ', '/opt/homebrew/bin/python3 ')) and 'train_teacher_student.py' in line for line in live)
    history=json.loads((OUT/'teacher_training_history.json').read_text())
    path=DATA/'models/teacher_student_selected.npz'
    receipt=json.loads((OUT/'teacher_data_receipt.json').read_text())
    sys.path.insert(0,str(OLD));import model
    base=model.read_weights()
    with np.load(path) as z:
        phase=str(z['selected_phase']);epoch=int(z['selected_epoch'])
        selected=next(r for r in history if r['phase']==phase and r['epoch']==epoch)
        keys=[k.removeprefix('encoder::') for k in z.files if k.startswith('encoder::')]
        changed={k:int(np.count_nonzero(z['encoder::'+k]!=base[k])) for k in keys}
        assert all(np.isfinite(z[k]).all() for k in z.files if np.issubdtype(z[k].dtype,np.number))
    assert selected['validation']['mean_within_anchor_Spearman']==max(r['validation']['mean_within_anchor_Spearman'] for r in history)
    frozen=json.loads((OUT/'frozen_teacher_checkpoint_receipt.json').read_text())
    assert sha(WORK/'train_teacher_student.py')==frozen['actual_running_training_source_sha256']
    assert sha(frozen['checkpoint'])==frozen['sha256']
    assert sha(OLD/'hyenadna/model.safetensors')==receipt['base_weights_sha256']
    assert sha(V3/'final_model_registry.json')==json.loads((OUT/'teacher_region_plan.json').read_text())['old_V3_registry_sha256']
    warmup=next(r for r in history if r['phase']==frozen['phase'] and r['epoch']==frozen['epoch'])
    result={
        'status':'Interrupted run; audited previously saved completed epoch only',
        'interruption':'Process no longer live. Log ends at native epoch2 batch600; epoch2 has no completed validation or eligible checkpoint.',
        'selected_phase':phase,'selected_epoch':epoch,
        'best_novel_anchor_Spearman':selected['validation']['mean_within_anchor_Spearman'],
        'selected_validation':selected['validation'],'warmup_validation':warmup['validation'],
        'selected_changed_native_scalars':sum(changed.values()),
        'selected_changed_native_scalars_by_parameter':{k:v for k,v in changed.items() if v},
        'checkpoint':str(path),'checkpoint_sha256':sha(path),
        'training_rows':sum(r['training_rows'] for r in receipt['datasets']),
        'new_exact_teacher_pairs':receipt['new_exact_teacher_rows'],
        'completed_native_epochs':sum(r['phase']=='native' for r in history),
        'reserved_teacher_test_used':False,'old_V3_registry_unchanged':True,
        'scope':'Imitates permitted AlphaGenome AVI predictions; no measured biological improvement or patient disease prediction established.',
        'next':'Matched measured-laboratory development comparison; all validation used here is reused development evidence.'}
    write_json(OUT/'teacher_training_status.json',result)
    print(json.dumps({k:result[k] for k in ['status','selected_phase','selected_epoch','selected_changed_native_scalars','checkpoint_sha256']}))
if __name__=='__main__':main()
