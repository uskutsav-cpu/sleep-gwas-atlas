from pathlib import Path
import csv,json,hashlib
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score

ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OLD=ROOT/'work/sleep_frailty_ml';OUT=ROOT/'outputs/sleep_frailty_ai_v2';MODELS=OUT/'models'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def metrics(y,p):
    y=np.asarray(y,dtype='int32');p=np.asarray(p,dtype='float64')
    assert y.shape==p.shape and np.isfinite(p).all()
    return {'n':len(y),'detected_calls':int(y.sum()),'prevalence':float(y.mean()),'auroc':float(roc_auc_score(y,p)),'average_precision':float(average_precision_score(y,p)),'brier_uncalibrated':float(np.mean((p-y)**2))}
def load():
    # Eager materialization avoids repeatedly decompressing entire arrays per batch.
    with np.load(OLD/'prepared.npz') as z:data={k:z[k] for k in z.files}
    with np.load(WORK/'sequence_features.npz') as z:seq={k:z[k] for k in z.files}
    audit=json.loads((OUT/'avi_retrieval_audit.json').read_text());assert not audit['partial_run'] and not audit['failed_groups']
    assert audit['output_sha256']==sha(OUT/'avi_exact_allele_scores.csv')
    manifest=list(csv.DictReader((ROOT/'outputs/sleep_frailty_ai_v1/split_manifest.csv').open()));assert len(manifest)==len(data['label'])
    avi={r['variant_GRCh38']:r for r in csv.DictReader((OUT/'avi_exact_allele_scores.csv').open())};aa=[]
    for row,y,c in zip(manifest,data['label'],data['chromosome']):
        assert int(row['label'])==int(y) and int(row['chromosome'])==int(c)
        chrom,pos,ref,alt=row['variant'].split('-');key=f'chr{chrom}:{pos}:{ref}:{alt}';r=avi.get(key)
        assert r is not None,('query never completed',key)
        aa.append([float(r['avi_raw']),float(r['avi_phred']),1.] if r['status']=='exact_allele_match' else [np.nan,np.nan,0.])
    aa=np.array(aa,dtype='float32');chrom=data['chromosome'];tr=np.flatnonzero(chrom<=16);va=np.flatnonzero((chrom>=17)&(chrom<=18));te=np.flatnonzero(chrom>=19)
    return data,seq,aa,manifest,tr,va,te
