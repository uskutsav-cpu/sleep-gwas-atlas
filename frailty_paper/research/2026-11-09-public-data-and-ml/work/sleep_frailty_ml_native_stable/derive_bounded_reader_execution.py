"""Create a separate execution revision, preserving both earlier attempts."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OLD=Path(__file__).parent
DEST=ROOT/'work/sleep_frailty_ml_bounded'
OUT=ROOT/'outputs/sleep_frailty_ml_bounded'


def main():
    DEST.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    assert not (OUT/'source_derivation_receipt.json').exists()
    common='''import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_native_stable'))
from stable_common import *
STABLE=DEST
DEST=ROOT/'outputs/sleep_frailty_ml_bounded'
original_assert_versions=assert_versions_stable

def assert_versions_stable():
    original_assert_versions()
    receipt=json.loads((DEST/'source_derivation_receipt.json').read_text())
    for path,h in receipt['preserved_previous_sources_and_attempts'].items():assert sha(ROOT/path)==h,path

def get_stable_plan():
    plan=json.loads((STABLE/'finite_stable_native_plan.json').read_text())
    plan['bounded_native_reader_execution']={
       'prior_CPU_plan_sha256':sha(STABLE/'finite_stable_native_plan.json'),
       'CPU_stop_receipt_sha256':sha(STABLE/'incomplete_CPU_preflight_stop.json'),
       'native_microbatch':32,
       'unused_rows':'Fill a short native-reader microbatch only with exact copies of its last real source construct. Remove those extra outputs before feature construction, head evaluation and losses. No DNA bases, flanks, measured targets or supervised observations are invented.',
       'equivalence':'Reader BN stays in eval mode and has no dropout; no operation mixes observations. Head is called once on the original complete actual observation batch. Discarded copies contribute no gradient.',
       'resource_policy':'MPS high/low watermark 0.4/0.3; real native reader input always 64 allele sequences by 4 bases by actual DNA length. Full 1000-update sustained preflight, all real length/tail/pair shapes, final complete validation pass.',
       'scientific_hyperparameters_loss_data_and_epochs_unchanged':True}
    plan['resource_revision']['device']='mps'
    path=DEST/'finite_stable_native_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    return plan
'''
    (DEST/'bounded_common.py').write_text(common)
    source=(OLD/'stable_network.py').read_text().replace('from stable_common import *','from bounded_common import *')
    begin=source.index('            selected=np.flatnonzero(packed.length[unique]==length);which=unique[selected];tokens=')
    end=source.index('        take=torch.from_numpy(inverse)',begin)
    replacement='''            selected=np.flatnonzero(packed.length[unique]==length);which=unique[selected]
            all_summaries=[];all_channels=[]
            for offset in range(0,len(which),32):
                real=which[offset:offset+32];nreal=len(real)
                # Exact real-construct copies bound the batch dimension only.
                # No extension or padding of the biological sequence is used.
                used=np.r_[real,np.repeat(real[-1],32-nreal)]
                tokens=packed.batch(used,device,dtype)
                if swap:tokens=torch.cat([tokens[32:],tokens[:32]])
                with torch.no_grad():prefix=self.reader.main[:3](self.reader.stem(tokens))
                activations=self.reader.main[3](prefix);g=activations.mean(-1)
                pos=torch.from_numpy(np.r_[packed.position[used],packed.position[used]].astype('int64')).to(device)
                local=local_pool(activations,pos);p=self.reader.head(self.reader.mapper(activations).mean(-1)).squeeze(-1)
                summary,blocks=summarize(g,local,p,32)
                all_summaries.append(summary[:nreal]);all_channels.append(((blocks@self.projection)*self.projection_multiplier)[:nreal])
            summary=torch.cat(all_summaries);projected=torch.cat(all_channels)
            index=torch.from_numpy(selected).to(device);summaries=summaries.index_copy(0,index,summary);channels=channels.index_copy(0,index,projected)
'''
    source=source[:begin]+replacement+source[end:]
    (DEST/'bounded_network.py').write_text(source)
    for filename in ['train_stable_development.py','run_stable_comparison.py','audit_fitted_stable.py','audit_stable_runtime.py']:
        s=(OLD/filename).read_text().replace('from stable_common import *','from bounded_common import *').replace('from stable_network import','from bounded_network import')
        if filename=='train_stable_development.py':s=s.replace("with_name('stable_network.py')","with_name('bounded_network.py')")
        if filename=='audit_stable_runtime.py':
            marker='    del old,new,a,b,p,p2\n'
            assert marker in s
            check='''    # Independent original full-batch versus bounded-reader gradients on
    # true neuron measurements. Both use double precision, fixed BN and no
    # head dropout. Extra unscored copies must have exactly zero influence.
    gradient_ids=np.r_[tr[neural[tr]][:3],tr[neural[tr]][:2]]
    gb={k:torch.from_numpy(d[k][gradient_ids]).double() for k in ['effect','signed_effect','activity_reference','activity_alternate']}
    gb['activity_center']=torch.from_numpy(centers[gradient_ids]).double()
    gw=torch.from_numpy(w[gradient_ids]).double();gq=torch.from_numpy(quality[gradient_ids]).double()
    old.zero_grad();new.zero_grad()
    original_loss(old(gradient_ids,d,packed,bank,x,extra,support,ls),gb,gw,gq,len(gradient_ids),0,0).backward()
    activity_loss(new(gradient_ids,d,packed,bank,x,extra,support,ls),gb,gw,gq,len(gradient_ids),0,0).backward()
    gradient_checks=[]
    for name,parameter in old.named_parameters():
        if not parameter.requires_grad:continue
        other=dict(new.named_parameters())[name]
        ga=parameter.grad if parameter.grad is not None else torch.zeros_like(parameter)
        gz=other.grad if other.grad is not None else torch.zeros_like(other)
        error=float((ga-gz).abs().max());scale=max(float(ga.abs().max()),float(gz.abs().max()),1e-8)
        assert error<1e-8+1e-7*scale,(name,error,scale)
        gradient_checks.append({'parameter':name,'maximum_absolute_error':error,'gradient_scale':scale})
    write_json(DEST/'bounded_calculation_equivalence.json',{'passed':True,'loss_checks':loss_checks,
        'duplicate_batch_forward_max_error':forward_error,'original_vs_bounded_gradient_checks':gradient_checks,
        'unscored_real_construct_copies_do_not_change_supervision':True,'new_test_outcomes_read':False})
    print('BOUNDED NATIVE CALCULATION EQUIVALENCE PASS',flush=True)
'''
            s=s.replace(marker,check+marker)
            s=s.replace("if (step+1)%50==0:","if step==0 or (step+1)%10==0:")
            s=s.replace("'code_sha256':{p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}}","'code_sha256':{p.name:sha(p) for p in Path(__file__).parent.glob('*.py')},'original_vs_bounded_gradient_checks':gradient_checks}")
        compile(s,str(DEST/filename),'exec');(DEST/filename).write_text(s)
    protected=[ROOT/'outputs/sleep_frailty_ml_native_stable/finite_stable_native_plan.json',
       ROOT/'outputs/sleep_frailty_ml_native_stable/incomplete_CPU_preflight_stop.json',
       ROOT/'outputs/sleep_frailty_ml_native_stable/mps_preflight_failure.json',
       OLD/'stable_network.py',OLD/'stable_loss.py',OLD/'stable_batches.py',OLD/'audit_stable_runtime.py']
    receipt={'preserved_previous_sources_and_attempts':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected},
       'new_code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DEST.glob('*.py')},
       'scientific_fits_started':False,'new_test_outcomes_used':False,'purpose':'Bound native-reader graph shapes while preserving all original observations, actual DNA lengths, complete main-batch objective, dropout location and training rules.'}
    (OUT/'source_derivation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('BOUNDED READER EXECUTION SOURCES DERIVED',flush=True)


if __name__=='__main__':main()
