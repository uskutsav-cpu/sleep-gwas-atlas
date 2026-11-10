"""Primal/dual bounds for the declared smooth absolute-error readout.

Iteration limits remain reported. A dual-feasible point gives an independent
lower bound for this fixed convex training objective, not an accuracy claim.
"""
from representation_common import *
from features import bank_features,design
from measurement_inputs import load_data
from training_utilities import transformations


def main():
    assert_versions_all();plan=get_plan();d,pairs,bank,x,extra=load_data()
    magnitude,_=bank_features(extra,bank,x,np.unique(d['bank_index']));results={}
    for source in plan['roles']:
        path=DEST/(source+'_diagnostic.json')
        if not path.exists(): continue
        record=json.loads(path.read_text());target=DEST/(source+'_readouts.npz');assert sha(target)==record['readouts_sha256']
        with np.load(target) as z: beta=z['magnitude_coefficient']
        train=d['split']==0
        if source!='none':
            sid=list(d['study_names']).index(source);groups=np.unique(d['group'][d['study']==sid]);train&=(d['study']!=sid)&~np.isin(d['group'],groups)
        ids=np.flatnonzero(train);assert hashlib.sha256(ids.tobytes()).hexdigest()==record['training_index_sha256']
        med,supported,weights,centers,quality,context=transformations(d,train,True)
        support=np.array(record['supported_lines']);line_med=np.array(record['line_medians']);nl=len(line_med);npen=len(beta)-nl
        # Replay the trainer's declared denominator exactly, including the
        # float32 weight reduction; the dual bound works for any positive a.
        a=weights[ids].astype('float64')/float(weights[ids].sum());lines=d['line'][ids];assert support[lines].all()
        base=np.where(support[lines],line_med[lines],med[d['cell'][ids]])
        residual=np.empty(len(ids))
        def blocks():
            for start in range(0,len(ids),8192):
                ix=ids[start:start+8192]
                yield start,design(magnitude,d['bank_index'][ix],d['line'][ix],support,record['magnitude_normalization'])
        for start,z in blocks():residual[start:start+len(z)]=base[start:start+len(z)]+z@beta-d['effect'][ids[start:start+len(z)]]
        root=np.sqrt(residual**2+.01**2);primal=float(a@root+.05*(beta[:npen]@beta[:npen]))
        assert abs(primal-record['magnitude_optimization']['final_objective'])<1e-9,(source,primal,record['magnitude_optimization']['final_objective'])
        u=residual/root
        for line in np.unique(lines):
            ix=lines==line;u[ix]-=float(a[ix]@u[ix]/a[ix].sum())
            u[ix]/=max(1.,float(abs(u[ix]).max()))
        feasible=max(abs(float(a[lines==line]@u[lines==line])) for line in np.unique(lines))
        assert feasible<1e-12 and abs(u).max()<=1+1e-12
        penalized_gradient=np.zeros(npen)
        for start,z in blocks():penalized_gradient+=z[:,:npen].T@(a[start:start+len(z)]*u[start:start+len(z)])
        t=base-d['effect'][ids]
        dual=float(a@(u*t+.01*np.sqrt(np.maximum(1-u*u,0)))-(penalized_gradient@penalized_gradient)/.2)
        gap=primal-dual;assert gap>=-1e-9
        results[source]={'primal':primal,'dual_lower_bound':dual,'certified_training_objective_gap_upper_bound':max(gap,0.),
            'fraction_of_primal':max(gap,0.)/primal,'maximum_dual_line_constraint_error':feasible,
            'optimizer_reported_success':record['magnitude_optimization']['success'],
            'predictive_validity_established':False}
    write_json(DEST/'convex_optimality_audit.json',{'fixed_objective_only':True,'completed_roles':results,
       'derivation':'For f(r)=sqrt(r^2+epsilon^2), f*(u)=-epsilon*sqrt(1-u^2), |u|<=1. Enforce weighted zero dual mean in each unpenalized line-offset column. Minimize the penalized slope term analytically to obtain the lower bound.',
       'regularization_lambda':.1,'smoothing_epsilon':.01,'bounds_do_not_establish_minimum_true_MAE_or_external_accuracy':True,
       'new_reserved_outcomes_read':False,'additional_model_fitting_or_penalty_search':False})
    assert_versions_all();print('CONVEX OBJECTIVE BOUNDS VERIFIED',len(results),max(r['certified_training_objective_gap_upper_bound'] for r in results.values()),flush=True)


if __name__=='__main__':main()
