"""One explicitly documented development-only loss correction before final testing."""
import shutil,subprocess,sys,argparse
from train import *
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-for-pid',type=int);args=parser.parse_args()
    if args.wait_for_pid:
        while subprocess.run(['ps','-p',str(args.wait_for_pid)],capture_output=True).returncode==0:time.sleep(10)
    assert not (OUT/'reserved_evaluation_started.json').exists() and not (OUT/'final_model_and_prediction_lock.json').exists()
    for name in ['development_selection.json','development_family_predictions.npz']:
        path=OUT/name
        if path.exists():
            target=path.with_name(path.stem+'_round1'+path.suffix);assert not target.exists();path.replace(target)
    plan=freeze_plan();initial=OUT/'model_search_plan_initial.json'
    if not initial.exists():
        shutil.copyfile(OUT/'model_search_plan.json',initial)
        receipt={'time_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'reason':'First-round measured-magnitude models can improve ranking but have negative macro raw-MAE skill versus a cell median. The log-Huber training objective does not directly minimize the prioritized reported absolute-size error. Add one finite raw-MAE correction round before opening any final-study model-performance result.','round1_records':{p.name:sha(p) for p in OUT.glob('development_*_record.json')},'round1_selection_sha256':sha(OUT/'development_selection_round1.json') if (OUT/'development_selection_round1.json').exists() else None,'interrupted_old_final_refits':'One original-recipe final checkpoint was already completed; another had one epoch when the parent was stopped. They used no final-study outcomes. Retain completed record and source curve as audit artifacts; unchanged recipes may reuse completed fits, partial fits restart from their fixed seed.','initial_plan_sha256':sha(initial),'additional_neural_fits':9,'architectures_seeds_rows_and_epoch_caps_changed':False,'final_study_model_performance_seen':False,'primary_selection_policy':'Preserve initial composite for comparisons. Also report whether any learned model beats the simple median on primary magnitude error; ranking gain alone is not evidence of stronger size prediction.'};write_json(OUT/'development_loss_amendment.json',receipt)
        plan['recipes'].update({'raw_MAE':{'pair_weight':0.,'no_cell':False,'raw_mae':True},'raw_MAE_pairs':{'pair_weight':.15,'no_cell':False,'raw_mae':True},'raw_MAE_no_cell_pairs':{'pair_weight':.15,'no_cell':True,'raw_mae':True}});plan['development_loss_amendment_sha256']=sha(OUT/'development_loss_amendment.json');write_json(OUT/'model_search_plan.json',plan)
    else:assert 'raw_MAE' in plan['recipes']
    # Numeric differentiation of the newly amended actual-magnitude loss.
    rng=np.random.default_rng(8801);c=rng.normal(0,.3,(12,560)).astype('float32');v=np.abs(rng.normal(0,.25,(12,623))).astype('float32');p={k:mx.array(a) for k,a in initialize(560,623,42).items()};u=mx.array(np.arange(12,dtype='int32')%2);target=mx.array(np.linspace(.2,1.2,12,dtype='float32'));scales=mx.array([.3,.2]);w=mx.ones(12);masks=[mx.ones((12,128)),mx.ones((12,64)),mx.ones((12,128)),mx.ones((12,64))]
    def objective(pp):return loss(pp,mx.array(c),mx.array(v),u,target,w,masks,4,2,2,.15,True,scales)
    value,g=mx.value_and_grad(objective)(p);direction=rng.normal(size=(64,2)).astype('float32');direction/=np.linalg.norm(direction);eps=.001;plus=dict(p);minus=dict(p);plus['effect_w']=p['effect_w']+eps*mx.array(direction);minus['effect_w']=p['effect_w']-eps*mx.array(direction);numeric=(float(objective(plus))-float(objective(minus)))/(2*eps);analytic=float(np.sum(np.array(g['effect_w'])*direction));assert abs(numeric-analytic)<max(.001,.05*abs(analytic)),(numeric,analytic);write_json(OUT/'raw_MAE_gradient_audit.json',{'passed':True,'numeric_gradient':numeric,'backprop_gradient':analytic,'continuous_raw_effect_error_loss':True,'final_study_performance_used':False})
    d=load();fit=np.flatnonzero(d['split']==0);va=np.flatnonzero(d['split']==1);mx.metal.set_memory_limit(1024*1024**2);mx.metal.set_cache_limit(64*1024**2)
    for recipe in ['raw_MAE','raw_MAE_pairs','raw_MAE_no_cell_pairs']:
        for seed in plan['seeds']:fit_nn(d,fit,va,recipe,seed)
    print('FINITE RAW-MAGNITUDE LOSS CORRECTION COMPLETE',flush=True)
if __name__=='__main__':main()
