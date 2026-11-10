"""Refit fixed screening-selected classical baselines on a development fold."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,pickle
from sklearn.base import clone
from common import *
assert_no_reserved_evaluation()
p=argparse.ArgumentParser();p.add_argument('--validation-chromosomes',required=True);a=p.parse_args();val=[int(c) for c in a.validation_chromosomes.split(',')]
plan=json.loads((OUT/'development_robustness_locked_config.json').read_text());selected=plan['classical_records'];data,features,manifest,_,_,_=load();tr,va=development_split(data,val);sid=data['study'];y=data['label']
path=OUT/'avi_exact_allele_scores.csv';audit=json.loads((OUT/'avi_retrieval_audit.json').read_text());assert not audit['failed_groups'] and audit['output_sha256']==sha(path)
lookup={r['variant_GRCh38']:r for r in csv.DictReader(path.open())};avi=np.array([[float(lookup[r['variant_GRCh38']]['avi_raw']),float(lookup[r['variant_GRCh38']]['avi_phred']),1.] if lookup[r['variant_GRCh38']]['status']=='exact_allele_match' else [np.nan,np.nan,0.] for r in manifest],dtype='float32')
with np.load(WORK/'reference_chromatin_features_local.npz') as z:chrom=z['features']
mat={'dna_pattern':features['dna_pattern'],'hyena_original':features['hyena_original'],'avi':avi,'dna_avi':np.c_[features['dna_pattern'],avi],'chromatin':chrom,'dna_chromatin':np.c_[features['dna_pattern'],chrom],'dna_chromatin_avi':np.c_[features['dna_pattern'],chrom,avi]}
prefix='development_baselines_val'+'-'.join(map(str,val));results=[];predictions={}
for source in selected:
    record=dict(source);name=prefix+'_'+source['name'];record['name']=name;record['screening_source_record']=source['name']
    if source['family']=='direct_avi':
        med=float(np.nanmedian(avi[tr,0]));v=np.where(np.isfinite(avi[va,0]),avi[va,0],med);pred=1/(1+np.exp(-np.clip(v.astype('float64'),-30,30)));record['config']={**source['config'],'imputation_raw_training_median':med};record['artifact']=None
    else:
        with (OUT/source['artifact']).open('rb') as f:old=pickle.load(f)
        models={};pred=np.empty(len(va))
        for s in range(2):
            tt=tr[sid[tr]==s];vv=np.flatnonzero(sid[va]==s);m=clone(old[s]);weight=source['config'].get('weighting','none');cw=None if weight=='none' else {0:1.,1:float(np.sqrt((len(tt)-y[tt].sum())/y[tt].sum()))}
            # Refitting a cloned pipeline recomputes imputation/scaling on this training fold.
            if source['family']=='study_logistic':
                step=list(m.named_steps)[-1];m.set_params(**{step+'__class_weight':cw});m.fit(mat[source['features']][tt],y[tt])
            else:
                weight=None if cw is None else np.where(y[tt]>0,cw[1],1.);kw={('histgradientboostingclassifier__sample_weight' if hasattr(m,'named_steps') else 'sample_weight'):weight};m.fit(mat[source['features']][tt],y[tt],**kw)
            pred[vv]=m.predict_proba(mat[source['features']][va[vv]])[:,1];models[s]=m
        target=MODELS/(name+'.pkl')
        with target.open('wb') as f:pickle.dump(models,f)
        record['artifact']=str(target.relative_to(OUT))
    record['validation']=study_metrics(y[va],pred,sid[va]);results.append(record);predictions[name]=pred;print('DEVELOPMENT BASELINE',name,record['validation']['macro_AP_skill'],flush=True)
(OUT/(prefix+'_search_results.json')).write_text(json.dumps(results,indent=2)+'\n');np.savez_compressed(OUT/(prefix+'_validation_predictions.npz'),row_index=va,label=y[va],study=sid[va],**predictions)
(OUT/(prefix+'_selection.json')).write_text(json.dumps({'status':'Fixed fold baselines complete','test_evaluated':False,'training_chromosomes':sorted(set(data['chromosome'][tr].tolist())),'validation_chromosomes':val},indent=2)+'\n')
