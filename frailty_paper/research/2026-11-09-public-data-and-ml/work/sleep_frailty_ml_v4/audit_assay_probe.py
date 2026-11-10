"""Replay the saved development ensembles and independently check metrics."""
import sys,subprocess
from common import *
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v3'))
from blocked_metrics import paired_block_bootstrap

def main():
    report=json.loads((OUT/'assay_transfer_results.json').read_text())
    status=json.loads((OUT/'teacher_training_status.json').read_text());teacher=Path(status['checkpoint'])
    assert sha(teacher)==report['teacher_checkpoint_sha256']
    with np.load(teacher) as z:fm=z['feature_mean'];fs=z['feature_sd'];hw=z['head.w1'];hb=z['head.b1']
    with np.load(OUT/'assay_transfer_validation_predictions.npz') as z:pred={k:z[k] for k in z.files}
    lookup={r['variant_GRCh38']:r for r in csv.DictReader((V3/'avi_exact_allele_scores.csv').open())}
    random=np.random.default_rng(41207).normal(0,1/np.sqrt(512),(512,64)).astype('float32')
    arrays={'original_control':[],'AVI_taught':[]};blocks=[];row_ids=[];labels=[];active=[];offset=0
    for stem in ['salomon','siraj']:
        sp=ROOT/'work/sleep_frailty_ml_v3'/(stem+'_sequences_local.npz')
        with np.load(sp) as z:ch=z['chromosome'];y=z['label']
        keep=np.flatnonzero(ch>=4);va=np.flatnonzero((ch[keep]>=4)&(ch[keep]<=5));ids=keep[va]
        row_ids.extend((offset+va).tolist());offset+=len(keep);labels.extend(y[ids].tolist())
        with np.load(sp.with_name(stem+'_features_local.npz')) as z:old=z['hyena_original'][ids]
        with np.load(DATA/(stem+'_taught_assay_features.npz')) as z:
            assert np.array_equal(z['row_index'],keep);assert str(z['signature'])==sha(teacher)+'_'+sha(sp)
            new=z['features'][va]
        rows=list(csv.DictReader(sp.with_name(stem+'_manifest_local.csv').open()));avi=[]
        for i in ids:
            r=rows[i];q=lookup[r['variant_GRCh38']];avi.append([float(q['avi_raw']),float(q['avi_phred']),1.] if q['status']=='exact_allele_match' else [np.nan,np.nan,0.])
            blocks.append(str(r['chromosome'])+':'+str(int(r['position_hg38'])//1000000));active.append(int(r['active_screen']))
        arrays['original_control'].append(np.c_[old,np.tanh(np.clip((old-fm)/fs,-20,20)@random),np.array(avi,dtype='float32')])
        arrays['AVI_taught'].append(np.c_[new,np.tanh(np.clip((new-fm)/fs,-20,20)@hw+hb),np.array(avi,dtype='float32')])
    assert np.array_equal(row_ids,pred['row_index']) and np.array_equal(labels,pred['label']) and np.array_equal(active,pred['active'])
    replay=[]
    for name,parts in arrays.items():
        x=np.concatenate(parts);members=[]
        for r in report['results']:
            if r['variant']!=name or r['seed']=='ensemble':continue
            assert sha(r['checkpoint'])==r['checkpoint_sha256']
            with np.load(r['checkpoint']) as z:p={k:z[k] for k in z.files}
            mean=p['feature_mean'];sd=p['feature_sd'];xx=np.clip((np.where(np.isfinite(x),x,mean)-mean)/sd,-20,20).astype('float32')
            def gelu(a):return .5*a*(1+np.tanh(.7978845608028654*(a+.044715*a**3)))
            # Match the saved inference batch size while using NumPy rather
            # than the MLX training implementation.
            scores=[]
            for start in range(0,len(xx),512):
                h=gelu(xx[start:start+512]@p['w1']+p['b1']);h=gelu(h@p['w2']+p['b2']);logit=(h@p['class.w']+p['class.b'])[np.arange(len(h)),pred['study'][start:start+len(h)]]
                scores.extend((1/(1+np.exp(-np.clip(logit,-30,30)))).tolist())
            members.append(scores)
        delta=float(np.max(np.abs(np.mean(members,axis=0)-pred[name])));assert delta<2e-5,(name,delta)
        replay.append({'model':name,'maximum_probability_error':delta,'checkpoint_count':len(members)})
    checks=[];uncertainty={}
    for s,stem in enumerate(['salomon','siraj']):
        for subset in ['all','active']:
            ix=np.flatnonzero((pred['study']==s)&((pred['active']==1) if subset=='active' else True))
            path=OUT/f'R_probe_{stem}_{subset}_input.csv';dest=path.with_name(path.stem+'_metrics.csv')
            with path.open('w',newline='') as f:
                writer=csv.writer(f);writer.writerow(['row_id','label','original_control','AVI_taught'])
                writer.writerows(zip(ix,pred['label'][ix],pred['original_control'][ix],pred['AVI_taught'][ix]))
            subprocess.run(['Rscript',str(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R'),str(path),str(dest)],check=True,capture_output=True,text=True)
            for r in csv.DictReader(dest.open()):
                expected=[roc_auc_score(pred['label'][ix],pred[r['model']][ix]),average_precision_score(pred['label'][ix],pred[r['model']][ix])]
                actual=[float(r['auroc']),float(r['average_precision'])];error=float(np.max(np.abs(np.array(expected)-actual)));assert error<1e-12
                checks.append({'study':stem,'subset':subset,'model':r['model'],'max_metric_error':error})
            uncertainty[stem+'_'+subset]=paired_block_bootstrap(pred['label'][ix],{name:pred[name][ix] for name in arrays},np.array(blocks)[ix],repetitions=1000,seed=3191+s)
    assert sha(V3/'final_model_registry.json')==json.loads((OUT/'teacher_region_plan.json').read_text())['old_V3_registry_sha256']
    write_json(OUT/'assay_probe_replay_audit.json',{'passed':True,'ensemble_replay':replay,'independent_R_checks':checks,'old_V3_registry_unchanged':True,'scope':'All comparisons and intervals are reused development evidence; stopping used these outcomes. No independent test or clinical validation.'})
    write_json(OUT/'assay_probe_block_uncertainty.json',uncertainty)
    print('Independent ensemble replay, 16 R metric checks and development block intervals passed.')
if __name__=='__main__':main()
