"""Paired holdout-chromosome bootstrap; descriptive with only four clusters."""
from pathlib import Path
import csv,json
import numpy as np

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'outputs/sleep_frailty_ai_v1'
with (OUT/'test_predictions.csv').open() as f:rows=list(csv.DictReader(f))
y=np.array([int(r['source_detected_call']) for r in rows]);chrom=np.array([int(r['chromosome']) for r in rows])
models=['constant','kmer_logistic','frozen_hyenadna','finetuned_hyenadna','small_dna_cnn'];pred={m:np.array([float(r[m]) for r in rows]) for m in models}
def metrics(y,p):
    ix=np.argsort(-p,kind='stable');s=p[ix];yy=y[ix];ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1]
    tp=np.cumsum(yy)[ends];fp=ends+1-tp;rec=tp/y.sum();precision=tp/(ends+1)
    return np.array([np.trapezoid(np.r_[0,rec],np.r_[0,fp/(len(y)-y.sum())]),np.sum(np.diff(np.r_[0,rec])*precision)])

rng=np.random.default_rng(26061106);groups=np.unique(chrom);draws={m:[] for m in models};diffs={m:[] for m in ['kmer_logistic','frozen_hyenadna','small_dna_cnn']};skipped=0
for _ in range(2000):
    selected=rng.choice(groups,len(groups),replace=True);ix=np.concatenate([np.flatnonzero(chrom==c) for c in selected])
    if not 0<y[ix].sum()<len(ix):skipped+=1;continue
    mm={m:metrics(y[ix],p[ix]) for m,p in pred.items()}
    for m in models:draws[m].append(mm[m])
    for m in diffs:diffs[m].append(mm['finetuned_hyenadna']-mm[m])
def describe(a):
    a=np.array(a)
    return {'auroc_percentile_range':np.quantile(a[:,0],[.025,.975]).tolist(),'average_precision_percentile_range':np.quantile(a[:,1],[.025,.975]).tolist()}
result={'seed':26061106,'requested_draws':2000,'skipped':skipped,'groups':groups.tolist(),'model_ranges':{m:describe(a) for m,a in draws.items()},'paired_finetuned_minus_comparator_ranges':{m:describe(a) for m,a in diffs.items()},'limitations':'Only four chromosome clusters. Exchangeable independent chromosomes are a rough assumption, not established. Percentile ranges are a descriptive sensitivity, not precise guarantees; training, architecture/split selection, assay-study and pretraining uncertainty are excluded.'}
(OUT/'test_uncertainty.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
