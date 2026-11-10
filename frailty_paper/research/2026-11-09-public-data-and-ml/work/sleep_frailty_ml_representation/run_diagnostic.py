"""Finite convex readout diagnostics, keeping measured magnitude primary."""
import time
from representation_common import *
from features import bank_features,training_normalization,design
from measurement_inputs import load_data
from training_utilities import transformations,detailed_metrics
from scipy.optimize import minimize
from scipy.stats import spearmanr


def main():
    assert_versions_all();plan=get_plan();start=time.monotonic();d,pairs,bank,x,extra=load_data();load_signed_sidecar(d)
    magnitude,signed=bank_features(extra,bank,x,np.unique(d['bank_index']));results={}
    for source in plan['roles']:
        path=DEST/(source+'_diagnostic.json')
        if path.exists(): results[source]=json.loads(path.read_text());continue
        print('START FIXED READOUT DIAGNOSTIC',source,flush=True);train=d['split']==0
        if source=='none':test=d['split']==1
        else:
            sid=list(d['study_names']).index(source);groups=np.unique(d['group'][d['study']==sid])
            train&=(d['study']!=sid)&~np.isin(d['group'],groups);test=(d['study']==sid)&(d['split']!=2)
        tr,te=np.flatnonzero(train),np.flatnonzero(test);assert not set(d['group'][tr])&set(d['group'][te])
        med,support,weights,centers,quality,context=transformations(d,train,True)
        nl=len(d['cell_line_names']);line_support=np.zeros(nl,dtype=bool);line_median=np.zeros(nl)
        for line in range(nl):
            ix=train&(d['line']==line);line_support[line]=ix.any();parent=np.unique(d['cell'][d['line']==line]);assert len(parent)<=1
            line_median[line]=np.median(d['effect'][ix]) if ix.any() else med[parent[0]] if len(parent) else np.median(d['effect'][tr])
        baseline=np.where(line_support[d['line']],line_median[d['line']],med[d['cell']])
        norm=training_normalization(magnitude,np.unique(d['bank_index'][tr]));p=140+6*nl+nl
        penalty=np.r_[np.full(140+6*nl,.1),np.zeros(nl)];total_weight=float(weights[tr].sum());history=[]
        def blocks(ids,matrix,settings,is_signed=False):
            for begin in range(0,len(ids),8192):
                ix=ids[begin:begin+8192]
                yield ix,design(matrix,d['bank_index'][ix],d['line'][ix],line_support,settings,is_signed)
        def objective(beta):
            value=0.;gradient=np.zeros(p)
            for ix,z in blocks(tr,magnitude,norm):
                residual=baseline[ix]+z@beta-d['effect'][ix];root=np.sqrt(residual*residual+.01**2);w=weights[ix]
                value+=float(w@root);gradient+=z.T@(w*residual/root)
            return value/total_weight+.5*float(penalty@(beta*beta)),gradient/total_weight+penalty*beta
        initial=objective(np.zeros(p))[0]
        def callback(beta):history.append(float(objective(beta)[0]))
        fit=minimize(objective,np.zeros(p),jac=True,method='L-BFGS-B',callback=callback,
                     options={'maxiter':80,'maxfun':120,'ftol':1e-10,'gtol':1e-6,'maxls':20})
        final,gradient=objective(fit.x);assert np.isfinite(fit.x).all() and final<=initial+1e-9
        prediction=np.empty(len(te));offset=0
        for ix,z in blocks(te,magnitude,norm):
            prediction[offset:offset+len(ix)]=np.maximum(0,baseline[ix]+z@fit.x);offset+=len(ix)
        # A distinct ridge readout tests signed-effect signal without replacing
        # the actual observed-magnitude endpoint by a latent target.
        signed_train=tr[np.isfinite(d['signed_effect'][tr])];signed_test=te[np.isfinite(d['signed_effect'][te])]
        signed_support=np.zeros(nl,dtype=bool);signed_support[np.unique(d['line'][signed_train])]=True
        snorm=training_normalization(signed,np.unique(d['bank_index'][signed_train]),True);ps=99+3*nl
        gram=np.zeros((ps,ps));rhs=np.zeros(ps);sw=0.
        for begin in range(0,len(signed_train),8192):
            ix=signed_train[begin:begin+8192];z=design(signed,d['bank_index'][ix],d['line'][ix],signed_support,snorm,True)
            w=(weights[ix]*quality[ix]).astype('float64');gram+=z.T@(w[:,None]*z);rhs+=z.T@(w*d['signed_effect'][ix]);sw+=float(w.sum())
        coefficient=np.linalg.solve(gram/sw+.1*np.eye(ps),rhs/sw)
        signed_prediction=np.empty(len(te));offset=0
        for begin in range(0,len(te),8192):
            ix=te[begin:begin+8192];z=design(signed,d['bank_index'][ix],d['line'][ix],signed_support,snorm,True)
            signed_prediction[offset:offset+len(ix)]=z@coefficient;offset+=len(ix)
        inverse=np.full(len(d['effect']),-1,dtype='int32');inverse[te]=np.arange(len(te));signed_metrics={}
        for c in np.unique(d['cell'][signed_test]):
            ix=signed_test[d['cell'][signed_test]==c];y=d['signed_effect'][ix];pred=signed_prediction[inverse[ix]];zero=float(abs(y).mean());mae=float(abs(y-pred).mean())
            signed_metrics[str(d['cell_names'][c])]={'n':len(ix),'MAE':mae,'zero_effect_MAE':zero,'skill_vs_zero':1-mae/zero,
                'Spearman':float(spearmanr(y,pred).statistic) if np.ptp(pred)>0 else None,
                'direction_accuracy_nonzero_effects':float((np.sign(y[y!=0])==np.sign(pred[y!=0])).mean()) if (y!=0).any() else None}
        outputs={'magnitude':prediction,'signed_effect':signed_prediction,'activity_reference':np.zeros(len(te)),'activity_alternate':np.zeros(len(te))}
        metrics=detailed_metrics(d,te,outputs,med,support,pairs);metrics.pop('activity_auxiliary',None);skills=[]
        for c in np.unique(d['cell'][te]):
            ix=d['cell'][te]==c;error=float(abs(d['effect'][te][ix]-prediction[ix]).mean());be=float(abs(d['effect'][te][ix]-baseline[te][ix]).mean());skill=1-error/be
            metrics['cells'][str(d['cell_names'][c])].update(line_baseline_MAE=be,MAE_skill_vs_line_baseline=skill)
            if support[c]:skills.append(skill)
        metrics['macro_MAE_skill_vs_line_baseline_supported_cells']=float(np.mean(skills)) if skills else None
        arrays=DEST/(source+'_readouts.npz');np.savez_compressed(arrays,observation_index=te,magnitude=prediction,signed_effect=signed_prediction,
            magnitude_coefficient=fit.x,signed_coefficient=coefficient,signed_gram=gram/sw,signed_rhs=rhs/sw)
        r={'role':source,'training_rows':len(tr),'prediction_rows':len(te),'original_signed_training_rows':len(signed_train),
           'training_index_sha256':hashlib.sha256(tr.tobytes()).hexdigest(),'prediction_index_sha256':hashlib.sha256(te.tobytes()).hexdigest(),
           'magnitude_normalization':norm,'signed_normalization':snorm,'line_medians':line_median.tolist(),'supported_lines':line_support.tolist(),'signed_supported_lines':signed_support.tolist(),
           'magnitude_optimization':{'success':bool(fit.success),'message':str(fit.message),'iterations':int(fit.nit),'function_evaluations':int(fit.nfev),
                'initial_objective':initial,'final_objective':final,'max_gradient':float(abs(gradient).max()),'callback_objective_history':history},
           'primary_magnitude_metrics':metrics,'signed_diagnostic':signed_metrics,'readouts_sha256':sha(arrays),
           'new_outside_study_used':False,'scientific_goal_complete':False}
        write_json(path,r);results[source]=r;assert_versions_all()
        print('DONE FIXED READOUT',source,'magnitude_skill',metrics['macro_MAE_skill_vs_line_baseline_supported_cells'],'converged',bool(fit.success),flush=True)
    mean=float(np.mean([results[s]['primary_magnitude_metrics']['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in plan['roles'][1:]]))
    write_json(DEST/'representation_diagnostic_assessment.json',{'all_eight_roles_complete':True,'mean_seven_held_source_magnitude_skill':mean,
       'per_source':{s:results[s]['primary_magnitude_metrics']['macro_MAE_skill_vs_line_baseline_supported_cells'] for s in plan['roles'][1:]},
       'all_magnitude_optimizations_converged':all(r['magnitude_optimization']['success'] for r in results.values()),
       'fixed_penalty_no_hyperparameter_search':True,'new_outside_validation_established':False,
       'failure_does_not_prove_no_nonlinear_or_richer_feature_signal':True,'scientific_goal_complete':False,'elapsed_seconds':time.monotonic()-start})
    assert_versions_all();print('REPRESENTATION DIAGNOSTIC COMPLETE',mean,flush=True)


if __name__=='__main__':main()
