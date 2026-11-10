"""Fixed full-AVI and physical-chromatin comparisons, validation only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','4')
import csv,json,pickle,numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from common import *
assert_no_reserved_evaluation()

data,features,manifest,tr,va,te=load();y=data['label'];study=data['study'];MODELS.mkdir(exist_ok=True)
path=OUT/'avi_exact_allele_scores.csv';audit=json.loads((OUT/'avi_retrieval_audit.json').read_text());assert not audit['failed_groups'] and audit['output_sha256']==sha(path) and audit['manifest_sha256']==sha(WORK/'avi_query_pairs.csv')
lookup={r['variant_GRCh38']:r for r in csv.DictReader(path.open())};avi=[]
for r in manifest:
    a=lookup[r['variant_GRCh38']];assert a['status']!='query_failed';avi.append([float(a['avi_raw']),float(a['avi_phred']),1] if a['status']=='exact_allele_match' else [np.nan,np.nan,0])
avi=np.asarray(avi,dtype='float32');cr=json.loads((OUT/'reference_chromatin_feature_receipt.json').read_text());cp=WORK/'reference_chromatin_features_local.npz';assert cr['feature_sha256']==sha(cp)
with np.load(cp) as z:chrom=z['features']
x={'avi':avi,'dna_avi':np.c_[features['dna_pattern'],avi],'chromatin':chrom,'dna_chromatin':np.c_[features['dna_pattern'],chrom],'dna_chromatin_avi':np.c_[features['dna_pattern'],chrom,avi]}
config={'features':{k:v.shape[1] for k,v in x.items()},'logistic_C':1.,'positive_weighting':['none','sqrt_negative_over_positive'],'boost_depths':[3,5],'boost_trees':300,'boost_lr':.05,'boost_l2':1.,'boost_min_leaf':40,'boost_max_bins':64,'preprocessing':'training-only median imputation, including explicit availability; training-only scaling for logistic','selection':'validation macro normalized AP skill','AVI_only':'One logistic comparator and direct raw-score ranking, not an unrestricted hyperparameter search','test_evaluated':False,'avi_sha256':sha(path),'reference_chromatin_sha256':sha(cp)}
(OUT/'full_avi_classical_locked_config.json').write_text(json.dumps(config,indent=2)+'\n');results=[];predictions={}
if (OUT/'full_avi_classical_results.json').is_file():
    results=json.loads((OUT/'full_avi_classical_results.json').read_text())
    with np.load(OUT/'full_avi_classical_validation_predictions.npz') as z:
        assert np.array_equal(z['row_index'],va) and np.array_equal(z['label'],y[va]);predictions={r['name']:z[r['name']] for r in results}
    for r in results:
        assert r['validation']==study_metrics(y[va],predictions[r['name']],study[va]);assert r['artifact'] is None or (OUT/r['artifact']).is_file()
completed=set(predictions)
def record(name,family,feature,cfg,pred,artifact):
    if name in completed:
        assert np.max(np.abs(pred-predictions[name]))<1e-12;return
    row={'name':name,'family':family,'features':feature,'config':cfg,'validation':study_metrics(y[va],pred,study[va]),'artifact':artifact};results.append(row);predictions[name]=pred;(OUT/'full_avi_classical_results.json').write_text(json.dumps(results,indent=2)+'\n');np.savez_compressed(OUT/'full_avi_classical_validation_predictions.npz',row_index=va,label=y[va],study=study[va],**predictions);print('FULL AVI VALIDATION',name,json.dumps(row['validation']),flush=True)
median=float(np.nanmedian(avi[tr,0]));raw=np.where(np.isfinite(avi[va,0]),avi[va,0],median);direct=1/(1+np.exp(-np.clip(raw.astype('float64'),-30,30)));record('full_avi_direct','direct_avi','avi',{'imputation_raw_training_median':median,'scope':'Monotonic ranking, not calibrated probability'},direct,None)
for feature in x:
    for weighting in (['none'] if feature=='avi' else ['none','sqrt']):
        name='full_'+feature+'_logistic_'+weighting;models={};pred=np.empty(len(va),dtype='float64')
        if name not in completed:
            for sid in range(2):
                tt=tr[study[tr]==sid];vv=np.flatnonzero(study[va]==sid);pw=float(np.sqrt((len(tt)-y[tt].sum())/y[tt].sum()));cw=None if weighting=='none' else {0:1.,1:pw};m=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),LogisticRegression(C=1.,class_weight=cw,max_iter=1500,random_state=42));m.fit(x[feature][tt],y[tt]);models[sid]=m;pred[vv]=m.predict_proba(x[feature][va[vv]])[:,1]
            target=MODELS/(name+'.pkl')
            with target.open('wb') as f:pickle.dump(models,f)
            record(name,'study_logistic',feature,{'C':1.,'weighting':weighting},pred,str(target.relative_to(OUT)))
        if feature=='avi':continue
        for depth in [3,5]:
            name=f'full_{feature}_boost_{weighting}_d{depth}';models={};pred=np.empty(len(va),dtype='float64')
            if name in completed:continue
            for sid in range(2):
                tt=tr[study[tr]==sid];vv=np.flatnonzero(study[va]==sid);pw=float(np.sqrt((len(tt)-y[tt].sum())/y[tt].sum()));weight=None if weighting=='none' else np.where(y[tt]>0,pw,1.);m=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),HistGradientBoostingClassifier(max_depth=depth,learning_rate=.05,max_iter=300,min_samples_leaf=40,l2_regularization=1.,max_bins=64,early_stopping=False,random_state=42));m.fit(x[feature][tt],y[tt],histgradientboostingclassifier__sample_weight=weight);models[sid]=m;pred[vv]=m.predict_proba(x[feature][va[vv]])[:,1]
            target=MODELS/(name+'.pkl')
            with target.open('wb') as f:pickle.dump(models,f)
            record(name,'study_hist_boosting',feature,{'depth':depth,'weighting':weighting,'trees':300,'learning_rate':.05,'l2':1.,'min_leaf':40,'max_bins':64},pred,str(target.relative_to(OUT)))
(OUT/'full_avi_classical_finished.json').write_text(json.dumps({'candidate_count':len(results),'test_evaluated':False,'status':'Fixed validation search complete'},indent=2)+'\n');print('FULL AVI CLASSICAL SEARCH COMPLETE',flush=True)
