import contextlib,io,json,runpy
from pathlib import Path
import numpy as np
with contextlib.redirect_stdout(io.StringIO()): d=runpy.run_path('work/neuronal_mpra_kmer_baseline.py')
pred=np.asarray(d['pred']);truth=np.asarray(d['truth']);assert len(pred)==161 and int(truth.sum())==32
old=json.loads(Path('outputs/stage6_neuronal_mpra_kmer_baseline.json').read_text())
def metrics(y,s):
 order=np.argsort(-s,kind='stable');y=y[order];s=s[order];n=int(y.sum());neg=len(y)-n;ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1];tp=np.cumsum(y)[ends];fp=ends+1-tp;rec=tp/n;prec=tp/(tp+fp);return float(np.trapezoid(np.r_[0,rec],np.r_[0,fp/neg])),float(np.sum(np.diff(np.r_[0,rec])*prec))
auc,ap=metrics(truth,pred);assert abs(auc-old['test_auroc'])<1e-12 and abs(ap-old['test_average_precision'])<1e-12
rng=np.random.default_rng(20261028);pos=np.flatnonzero(truth==1);neg=np.flatnonzero(truth==0);b=[]
for _ in range(5000):
 inds=np.concatenate([rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)]);b.append(metrics(truth[inds],pred[inds]))
b=np.array(b);out={'n_resamples':len(b),'seed':20261028,'point_auc':auc,'auc_percentile_95':np.quantile(b[:,0],[.025,.975]).tolist(),'point_ap':ap,'ap_percentile_95':np.quantile(b[:,1],[.025,.975]).tolist(),'method':'Stratified resampling of test variants within observed outcome labels; fixed model and fixed 32/161 prevalence','limitations':'Conditional sampling variability only. Variant resampling assumes exchangeability and does not account for chromosome or regional LD clustering, training uncertainty, chosen chromosome split, selection or assay-specific detection labels; two held-out chromosomes are insufficient for robust cluster inference. Not disease-score validation.'}
Path('outputs/stage6_neuronal_mpra_kmer_uncertainty.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
