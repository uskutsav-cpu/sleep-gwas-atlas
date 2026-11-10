#!/usr/bin/env python3
"""Evaluate externally generated, fixed predictions; never train models or output rows."""
import argparse, csv, json, math, pathlib
import numpy as np
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.optimize import minimize
from scipy.special import expit


def read_data(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        needed = {'outcome', 'baseline_probability', 'augmented_probability'}
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)) or not needed.issubset(reader.fieldnames):
            raise ValueError('Need unique outcome, baseline_probability and augmented_probability columns')
        result=[]
        for n,r in enumerate(reader,2):
            try: v=[float(r[x]) for x in ['outcome','baseline_probability','augmented_probability']]
            except (ValueError,TypeError): raise ValueError(f'Invalid or missing numeric value on row {n}')
            if not all(math.isfinite(x) for x in v) or v[0] not in (0.,1.) or not all(0 < p < 1 for p in v[1:]):
                raise ValueError(f'Invalid outcome or probability on row {n}; probabilities must be strictly between 0 and 1')
            result.append(v)
    a=np.asarray(result,dtype=float)
    if len(a)<2 or len(np.unique(a[:,0]))!=2: raise ValueError('Both outcome classes are required')
    return a


def metrics(a):
    y,p,q=a.T
    auc0=roc_auc_score(y,p); auc1=roc_auc_score(y,q)
    b0=brier_score_loss(y,p); b1=brier_score_loss(y,q)
    return {'baseline_auc':float(auc0),'augmented_auc':float(auc1), 'auc_difference_augmented_minus_baseline':float(auc1-auc0),
            'baseline_brier':float(b0),'augmented_brier':float(b1),'brier_difference_augmented_minus_baseline':float(b1-b0)}


def calibration(y,p):
    x=np.log(p/(1-p))
    # Calibration-in-the-large uses a fixed slope 1; slope/intercept fitted jointly.
    def loss(v, offset, slopefree):
        eta=offset+v[0]+(v[1]*x if slopefree else 0)
        return float(np.sum(np.logaddexp(0,eta)-y*eta))
    def gradient(v, offset, slopefree):
        eta=offset+v[0]+(v[1]*x if slopefree else 0)
        residual=expit(eta)-y
        return np.array([residual.sum(), (residual*x).sum()] if slopefree else [residual.sum()])
    a=minimize(loss,[0.],args=(x,False),jac=gradient,method='BFGS')
    b=minimize(loss,[0.,1.],args=(np.zeros_like(x),True),jac=gradient,method='BFGS')
    if not a.success or not b.success or max(abs(b.x))>20 or np.ptp(x)<1e-10:
        return {'status':'not reliably estimable','calibration_in_the_large':None,'intercept_joint':None,'slope':None}
    return {'status':'numerically converged; assess uncertainty and separation','calibration_in_the_large':float(a.x[0]),'intercept_joint':float(b.x[0]),'slope':float(b.x[1])}


def evaluate(a, iterations, seed):
    if iterations<200: raise ValueError('At least 200 bootstrap draws required')
    point=metrics(a); samples={k:[] for k in point}; rng=np.random.default_rng(seed); invalid=0
    # Ordinary paired iid resampling; draws with one outcome class cannot compute AUC.
    for _ in range(iterations):
        b=a[rng.integers(0,len(a),len(a))]
        if len(np.unique(b[:,0]))<2: invalid+=1; continue
        vals=metrics(b)
        for k,v in vals.items(): samples[k].append(v)
    if invalid/iterations>0.01: raise ValueError('Over 1% of draws have one class; sample too sparse for this iid bootstrap')
    ci={k:[float(v) for v in np.quantile(s,[.025,.975])] for k,s in samples.items()}
    return {'n':len(a),'events':int(a[:,0].sum()),'metrics':point,'paired_iid_bootstrap_percentile_95':ci,
            'bootstrap_draws':iterations,'invalid_one_class_draws':invalid,'seed':seed,
            'baseline_calibration':calibration(a[:,0],a[:,1]),'augmented_calibration':calibration(a[:,0],a[:,2]),
            'limitations':['Assumes independent, equally sampled complete records; invalid for unhandled clusters, survey weights or censoring.',
             'Fixed predictions must come from independently fitted and frozen models; code cannot verify data provenance or selection overlap.',
             'Bootstrap uncertainty conditions on fixed model fitting; it does not include training or tuning uncertainty.',
             'Calibration estimates are descriptive for the evaluation sample and may be unstable in small samples.']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--bootstrap',type=int,default=2000);p.add_argument('--seed',type=int,default=20261029)
    p.add_argument('--confirm-independent-fixed-predictions',action='store_true',required=True,help='Analyst asserts, not proves, model and target independence')
    args=p.parse_args(); out=pathlib.Path(args.output)
    if out.exists(): p.error('Output already exists')
    try:
        a=read_data(args.input); result=evaluate(a,args.bootstrap,args.seed)
        with out.open('x',encoding='utf-8') as f: json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    except (ValueError,OSError) as e: p.error(str(e))

if __name__=='__main__':main()
