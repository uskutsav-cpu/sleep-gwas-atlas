"""Separate recovery version; scientific operations and finite gates unchanged."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OLD=Path(__file__).parent
NEW=ROOT/'work/sleep_frailty_ml_resumable'
OUT=ROOT/'outputs/sleep_frailty_ml_resumable'


def once(s,a,b):
    assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)


def main():
    NEW.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    assert not (OUT/'source_derivation_receipt.json').exists()
    common='''import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_bounded'))
from bounded_common import *
BOUNDED=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_resumable'
original_assert_versions=assert_versions_stable
sys.path.insert(0,str(Path(__file__).resolve().parent))

def assert_versions_stable():
    original_assert_versions()
    manifest=json.loads((BOUNDED/'interrupted_run_manifest.json').read_text())
    for path,h in manifest['files_sha256'].items():assert sha(ROOT/path)==h,path
    receipt=json.loads((DEST/'source_derivation_receipt.json').read_text())
    for path,h in receipt['previous_code_sha256'].items():assert sha(ROOT/path)==h,path

def get_stable_plan():
    plan=json.loads((BOUNDED/'finite_stable_native_plan.json').read_text())
    plan['recovery_execution']={
        'prior_plan_sha256':sha(BOUNDED/'finite_stable_native_plan.json'),
        'interruption_manifest_sha256':sha(BOUNDED/'interrupted_run_manifest.json'),
        'checkpoint_every_updates':100,
        'checkpoint_contents':'Trainable parameters, Adam moments/steps, CPU and MPS RNG, NumPy sampler state, saved epoch losses and completed curves; atomic replacement; weights_only=True loading.',
        'resumption':'Reconstruct and verify the original batch/pair schedule from seed without optimizer steps for past batches. Restore state before the next real update. Preserve artifacts beyond the last committed checkpoint before replay.',
        'initial_recovery':'Original interrupted run lacks optimizer state. This version starts again from the same pinned initial weights, preserving the partial original; it does not initialize from its partial weights or select between repeats.',
        'scientific_labels_model_objective_epochs_seed_and_gates_unchanged':True,
        'new_neuronal_activity_inventory_used_for_fit':False}
    path=DEST/'finite_stable_native_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan
'''
    (NEW/'resumable_common.py').write_text(common)
    (NEW/'resumption.py').write_bytes((OLD/'resumption.py').read_bytes())
    s=(OLD/'train_stable_development.py').read_text().replace('from bounded_common import *','from resumable_common import *')
    s=s.replace('import argparse,time,shutil','import argparse,time,shutil\nfrom resumption import save_recovery,load_recovery')
    s=once(s,"parser.add_argument('--held-study',required=True);args=parser.parse_args()","parser.add_argument('--held-study',required=True);parser.add_argument('--resume',action='store_true');args=parser.parse_args()")
    begin=s.index("    proof=json.loads((DEST/'stable_runtime_preflight.json')")
    end=s.index('    assert shutil.disk_usage(ROOT).free',begin)
    s=s[:begin]+'''    proof=json.loads((BOUNDED/'stable_runtime_preflight.json').read_text());assert proof['passed']
    assert json.loads((DEST/'recovery_state_audit.json').read_text())['passed']
    for filename in ['train_stable_development.py','bounded_network.py','bounded_common.py']:
        assert sha(ROOT/'work/sleep_frailty_ml_bounded'/filename)==proof['code_sha256'][filename]
'''+s[end:]
    s=once(s,"    assert not (directory/(name+'_started.json')).exists(),'Preserve and inspect any incomplete prior fit before a restart.'", "    assert (directory/(name+'_started.json')).exists()==args.resume,'A fresh fit and a resumed fit require distinct explicit state.'")
    s=s.replace("Path(__file__).with_name('bounded_network.py')","ROOT/'work/sleep_frailty_ml_bounded/bounded_network.py'")
    s=once(s,"    write_json(directory/(name+'_started.json'),meta);opt=plan['optimization']",'''    identity={k:meta[k] for k in ['recipe','held_study','seed','epochs','training_index_sha256','prediction_index_sha256','initial_parameter_sha256','plan_sha256','code_sha256']}
    identity['recovery_helper_sha256']=sha(Path(__file__).with_name('resumption.py'))
    if args.resume:
        original=json.loads((directory/(name+'_started.json')).read_text())
        for key in identity:
            if key!='recovery_helper_sha256':assert original[key]==meta[key],key
    else:write_json(directory/(name+'_started.json'),meta)
    opt=plan['optimization']''')
    old="    rng=np.random.default_rng(seed);schedule=hashlib.sha256();updates=0;curve=[];start=time.monotonic();last_print=start"
    new='''    recovery=directory/(name+'_recovery.pt')
    rng=np.random.default_rng(seed);schedule=hashlib.sha256();updates=0;curve=[];start=time.monotonic();last_print=start
    restored=load_recovery(recovery,model,optimizer,device,identity) if args.resume else None
    saved_updates=int(restored['updates']) if restored else 0;updates=saved_updates
    saved_epoch=int(restored['epoch']) if restored else 0
    elapsed_before=float(restored['elapsed_seconds']) if restored else 0.
    curve=restored['curve'] if restored else []
    if restored:
        # Preserve any writes made after the last atomically committed state.
        beyond=[]
        for path in directory.glob(name+'_epoch*'):
            number=int(path.name.split('_epoch',1)[1].split('_',1)[0].split('.',1)[0])
            if number>len(curve):beyond.append(path)
        curvepath=directory/(name+'_curve.json')
        if curvepath.exists():beyond.append(curvepath)
        if beyond:
            archive=directory/('recovery_archive_'+str(time.time_ns()));archive.mkdir()
            manifest={p.name:sha(p) for p in beyond}
            for path in beyond:path.rename(archive/path.name)
            write_json(archive/'preserved_after_last_committed_state.json',manifest)
        history=directory/(name+'_resume_history.json')
        events=json.loads(history.read_text()) if history.exists() else []
        events.append({'PID':os.getpid(),'recovery_sha256':sha(recovery),'updates':saved_updates,'epoch':saved_epoch})
        write_json(history,events)
        if curve:write_json(curvepath,curve)
    ordinal=0;chain_verified=saved_updates==0
    if restored and saved_updates==0:
        assert restored['schedule_sha256']==schedule.hexdigest()
        assert json.loads(restored['NumPy_rng_json'])==rng.bit_generator.state
    epoch=1;objectives=[];anchors=[]
    def checkpoint_now(phase):
        save_recovery(recovery,model,optimizer,device,identity,{'updates':updates,'epoch':epoch,'phase':phase,
           'schedule_sha256':schedule.hexdigest(),'NumPy_rng_json':json.dumps(rng.bit_generator.state),
           'objectives':objectives,'anchors':anchors,'curve':curve,'elapsed_seconds':elapsed_before+time.monotonic()-start})
    if not restored:
        torch.manual_seed(seed+1)
        if device=='mps':torch.mps.manual_seed(seed+1)
        checkpoint_now('training')'''
    s=once(s,old,new)
    old='''        torch.manual_seed(seed+epoch)
        if device=='mps':torch.mps.manual_seed(seed+epoch)
        model.train();batches=epoch_batches(tr,lengths,rng,opt['main_batch']);objectives=[];anchors=[];completed_rows=0'''
    new='''        if not restored or epoch>saved_epoch:
            torch.manual_seed(seed+epoch)
            if device=='mps':torch.mps.manual_seed(seed+epoch)
        model.train();batches=epoch_batches(tr,lengths,rng,opt['main_batch']);completed_rows=0
        objectives=list(restored['objectives']) if restored and epoch==saved_epoch and restored['phase']=='training' else []
        anchors=list(restored['anchors']) if restored and epoch==saved_epoch and restored['phase']=='training' else []'''
    s=once(s,old,new)
    old="            completed_rows+=len(main_ids);schedule.update(ids.tobytes());optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls)"
    new='''            completed_rows+=len(main_ids);schedule.update(ids.tobytes());ordinal+=1
            if ordinal<=saved_updates:
                if ordinal==saved_updates:
                    assert epoch==saved_epoch
                    assert schedule.hexdigest()==restored['schedule_sha256']
                    assert rng.bit_generator.state==json.loads(restored['NumPy_rng_json'])
                    chain_verified=True
                continue
            assert chain_verified and updates==ordinal-1
            optimizer.zero_grad();prediction=model(ids,d,packed,bank,x,extra,support,ls)'''
    s=once(s,old,new)
    s=once(s,"            if device=='mps' and updates%100==0:torch.mps.empty_cache()", "            if device=='mps' and updates%100==0:torch.mps.empty_cache()\n            if updates%plan['recovery_execution']['checkpoint_every_updates']==0:checkpoint_now('training')")
    s=once(s,"        if device=='mps':torch.mps.synchronize();torch.mps.empty_cache()", "        if restored and (epoch<saved_epoch or (epoch==saved_epoch and restored['phase']=='epoch_complete')):continue\n        if device=='mps':torch.mps.synchronize();torch.mps.empty_cache()")
    s=s.replace('time.monotonic()-start', 'elapsed_before+time.monotonic()-start')
    s=s.replace("'elapsed_seconds':elapsed_before+elapsed_before+time.monotonic()-start", "'elapsed_seconds':elapsed_before+time.monotonic()-start")
    s=once(s,"        print('NATIVE EPOCH COMPLETE'", "        checkpoint_now('epoch_complete')\n        print('NATIVE EPOCH COMPLETE'")
    compile(s,str(NEW/'train_resumable_development.py'),'exec');(NEW/'train_resumable_development.py').write_text(s)
    runner=(OLD/'run_stable_comparison.py').read_text().replace('from bounded_common import *','from resumable_common import *').replace("with_name('train_stable_development.py')","with_name('train_resumable_development.py')")
    runner=once(runner,"    assert not (directory/(name+'_started.json')).exists(),('Incomplete existing fit requires process/checkpoint inspection; do not silently restart',name)",'''    started=(directory/(name+'_started.json')).exists()
    if started:assert (directory/(name+'_recovery.pt')).exists(),('Explicit recovery requires saved full state',name)''')
    runner=once(runner,"with (directory/(name+'.log')).open('w') as handle:","with (directory/(name+'.log')).open('a' if started else 'w') as handle:")
    runner=once(runner,"'--recipe',recipe,'--held-study',source],stdout=handle", "'--recipe',recipe,'--held-study',source]+(['--resume'] if started else []),stdout=handle")
    runner=runner.replace("assert json.loads((DEST/'stable_runtime_preflight.json').read_text())['passed']", "assert json.loads((BOUNDED/'stable_runtime_preflight.json').read_text())['passed'];assert json.loads((DEST/'recovery_state_audit.json').read_text())['passed']")
    compile(runner,str(NEW/'run_resumable_comparison.py'),'exec');(NEW/'run_resumable_comparison.py').write_text(runner)
    receipt={'previous_code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD/'train_stable_development.py',OLD/'bounded_network.py',OLD/'bounded_common.py',OLD/'run_stable_comparison.py']},
       'new_code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in NEW.glob('*.py')},
       'changes':'Only recovery state persistence, explicit resume, schedule replay and preservation of uncommitted artifacts. Same model, losses, actual DNA, labels, batches, three epochs, seeds, optimization and finite neural gate.',
       'scientific_fit_started':False,'new_test_outcomes_used':False}
    (OUT/'source_derivation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('RESUMABLE FINITE COMPARISON SOURCES DERIVED',flush=True)


if __name__=='__main__':main()
