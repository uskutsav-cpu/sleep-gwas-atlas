from pathlib import Path
import csv,json,hashlib
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score,r2_score
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v3';OLD=ROOT/'work/sleep_frailty_ml';MODELS=OUT/'models';STUDIES=['Salomon_NGN2','Siraj_OL42_SKNSH']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def assert_no_reserved_evaluation():
    assert not (OUT/'reserved_evaluation_started.json').exists(),'V3 weights/choices are frozen after the reserved evaluation starts. Preserve the result; do not refit this version.'
def metrics(y,p):
    y=np.asarray(y,dtype='int32');p=np.asarray(p,dtype='float64');assert y.shape==p.shape and np.isfinite(p).all()
    ap=float(average_precision_score(y,p));prev=float(y.mean())
    return {'n':len(y),'detected_calls':int(y.sum()),'prevalence':prev,'auroc':float(roc_auc_score(y,p)),'average_precision':ap,'AP_skill_over_prevalence':(ap-prev)/(1-prev),'brier_uncalibrated':float(np.mean((p-y)**2))}
def study_metrics(y,p,study):
    rr={name:metrics(y[study==i],p[study==i]) for i,name in enumerate(STUDIES)}
    return {'studies':rr,'macro_AP_skill':float(np.mean([r['AP_skill_over_prevalence'] for r in rr.values()]))}
def load():
    receipt=json.loads((OUT/'feature_receipt.json').read_text());arrays=[];features=[];manifests=[]
    for i,stem in enumerate(['salomon','siraj']):
        path=WORK/(stem+'_sequences_local.npz');fp=WORK/(stem+'_features_local.npz');r=next(r for r in receipt['datasets'] if r['study']==stem)
        assert r['input_sha256']==sha(path) and r['features_sha256']==sha(fp)
        with np.load(path) as z:d={k:z[k] for k in z.files}
        with np.load(fp) as z:f={k:z[k] for k in z.files}
        manifest_path=WORK/(stem+'_manifest_local.csv');rows=list(csv.DictReader(manifest_path.open()));assert len(rows)==len(d['label'])==len(f['hyena_original'])
        chromatin_receipt=OUT/'reference_chromatin_feature_receipt.json'
        if chromatin_receipt.exists():assert json.loads(chromatin_receipt.read_text())['manifest_sha256'][stem]==sha(manifest_path)
        for j,row in enumerate(rows):assert int(row['label'])==d['label'][j] and int(row['chromosome'])==d['chromosome'][j]
        arrays.append(d);features.append(f);manifests.extend(rows)
    data={k:np.concatenate([d[k] for d in arrays]) for k in ['label','chromosome','activity_mean','absolute_effect']};data['study']=np.concatenate([np.full(len(d['label']),i,dtype='int32') for i,d in enumerate(arrays)])
    ff={k:np.concatenate([f[k] for f in features]).astype('float32') for k in features[0]};c=data['chromosome'];tr=np.flatnonzero(c>=6);va=np.flatnonzero((c>=4)&(c<=5));te=np.flatnonzero(c<=3)
    assert set(tr).isdisjoint(va) and set(tr).isdisjoint(te) and set(va).isdisjoint(te)
    return data,ff,manifests,tr,va,te

def prior_avi(manifest):
    """Stable permitted-score snapshot for stage1 auxiliary supervision only."""
    known={};sources=[]
    for stem in ['avi','external_avi']:
        path=ROOT/'outputs/sleep_frailty_ai_v2'/(stem+'_exact_allele_scores.csv');receipt=json.loads((path.parent/(stem+'_retrieval_audit.json')).read_text());assert receipt['output_sha256']==sha(path)
        sources.append({'file':str(path),'sha256':sha(path)})
        for r in csv.DictReader(path.open()):
            if r['status']=='exact_allele_match':known[r['variant_GRCh38']]=float(r['avi_raw'])
    target=np.array([known.get(r['variant_GRCh38'],np.nan) for r in manifest],dtype='float32')
    return target,sources

def development_split(data,validation_chromosomes):
    """Development-only split; any extra folds require a sequence overlap audit."""
    val=sorted(set(validation_chromosomes));assert val and all(4<=c<=22 for c in val)
    c=data['chromosome'];tr=np.flatnonzero((c>=4)&~np.isin(c,val));va=np.flatnonzero(np.isin(c,val))
    if val!=[4,5]:
        path=OUT/'development_fold_sequence_audit.json';receipt=json.loads(path.read_text())
        for stem in ['salomon','siraj']:assert receipt['source_sha256'][stem]==sha(WORK/(stem+'_sequences_local.npz'))
        fold=next(f for f in receipt['folds'] if f['validation_chromosomes']==val)
        assert fold['validation_indices_sha256']==hashlib.sha256(va.tobytes()).hexdigest()
        tr=tr[~np.isin(tr,fold['excluded_training_indices'])]
        assert fold['training_indices_sha256']==hashlib.sha256(tr.tobytes()).hexdigest()
    assert len(tr) and len(va) and not set(tr)&set(va) and np.all(c[tr]>=4) and np.all(c[va]>=4)
    return tr,va
