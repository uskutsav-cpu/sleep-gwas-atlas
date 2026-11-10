"""Reload locally trained model records for batched research inference.

This module does not fit anything and does not select/evaluate a test set.
"""
import json,csv,pickle,numpy as np
from common import *
from contrast_features import ContrastFeatures
from contrast_network import forward as contrast_forward

class StoredFeatureInference:
    def __init__(self):
        self.features=ContrastFeatures();self.data=self.features.data;self.original=self.features.original;self.manifest=self.features.manifest;self.avi=None;self.cache={}
    def avi_features(self):
        if self.avi is None:
            path=OUT/'avi_exact_allele_scores.csv';audit=json.loads((OUT/'avi_retrieval_audit.json').read_text());assert not audit['failed_groups'] and audit['output_sha256']==sha(path) and audit['manifest_sha256']==sha(WORK/'avi_query_pairs.csv');lookup={r['variant_GRCh38']:r for r in csv.DictReader(path.open())};values=[]
            for r in self.manifest:
                row=lookup[r['variant_GRCh38']];assert row['status']!='query_failed';values.append([float(row['avi_raw']),float(row['avi_phred']),1.] if row['status']=='exact_allele_match' else [np.nan,np.nan,0.])
            self.avi=np.asarray(values,dtype='float32')
        return self.avi
    def motif(self,ix):
        f=self.features;m=np.empty((len(ix),3*f.M),dtype='float32');lo=ix<f.nsal;m[lo]=f.motifs[0][ix[lo]];m[~lo]=f.motifs[1][ix[~lo]-f.nsal];return m
    def matrix(self,kind,ix,avi_input=False):
        if kind=='dna_pattern':x=self.original['dna_pattern'][ix]
        elif kind=='hyena_original':x=self.original['hyena_original'][ix]
        elif kind=='avi':x=self.avi_features()[ix]
        elif kind=='dna_avi':x=np.c_[self.original['dna_pattern'][ix],self.avi_features()[ix]]
        elif kind=='chromatin':x=self.features.chromatin[ix]
        elif kind=='dna_chromatin':x=np.c_[self.original['dna_pattern'][ix],self.features.chromatin[ix]]
        elif kind=='dna_chromatin_avi':x=np.c_[self.original['dna_pattern'][ix],self.features.chromatin[ix],self.avi_features()[ix]]
        elif kind=='motif':x=self.motif(ix)
        elif kind=='motif_chromatin':x=np.c_[self.motif(ix),self.features.chromatin[ix]]
        elif kind=='hyena_chromatin':x=np.c_[self.original['hyena_original'][ix],self.features.chromatin[ix]]
        else:raise ValueError(kind)
        if avi_input:x=np.c_[x,self.avi_features()[ix]]
        return x.astype('float32')
    def predict(self,record,indices,registry):
        ix=np.asarray(indices);key=(record['name'],sha_indices(ix))
        if key in self.cache:return self.cache[key]
        family=record['family'];sid=self.data['study'][ix]
        if family in ['seed_ensemble','biological_seed_ensemble','context_variant_ensemble']:
            members=[self.predict(registry[n],ix,registry) for n in record['config']['members']];result={k:np.mean([v[k] for v in members],axis=0) for k in members[0]}
        elif family=='direct_avi':
            a=self.avi_features()[ix,0];a=np.where(np.isfinite(a),a,record['config']['imputation_raw_training_median']);result={'overall':sigmoid(a)}
        elif family in ['study_logistic','study_hist_boosting']:
            # Only our own local checkpoints; never load untrusted pickle files.
            with (OUT/record['artifact']).open('rb') as f:models=pickle.load(f)
            pred=np.empty(len(ix),dtype='float64')
            for s in range(2):
                use=np.flatnonzero(sid==s)
                if len(use):pred[use]=models[s].predict_proba(self.matrix(record['features'],ix[use]))[:,1]
            result={'overall':pred}
        elif family in ['frozen_multi_study_head','biological_feature_head']:
            z=np.load(OUT/record['artifact']);cfg=json.loads((OUT/record['config']['locked_config']).read_text()) if family=='biological_feature_head' else {};pred=[]
            for a in range(0,len(ix),512):
                jj=ix[a:a+512];s=self.data['study'][jj];x=self.matrix(record['features'],jj,cfg.get('AVI_as_input',False));x=np.clip((np.where(np.isfinite(x),x,z['feature_mean'])-z['feature_mean'])/z['feature_sd'],-20,20).astype('float32');h=gelu(x@z['shared.w1']+z['shared.b1']);h=gelu(h@z['shared.w2']+z['shared.b2']);logits=(h@z['class.w']+z['class.b'])[np.arange(len(jj)),s];pred.extend(sigmoid(logits))
            result={'overall':np.asarray(pred)}
        elif family=='context_variant_mixture':
            z=np.load(OUT/record['artifact']);params={k:z[k] for k in z.files if k.startswith(('context.','variant.','modulation.','gate.','conditional.','activity.','effect.','avi.'))};self.features.load_transform(z);values={k:[] for k in ['overall','conditional_scores','predicted_active','absolute_effect']}
            for a in range(0,len(ix),512):
                jj=ix[a:a+512];c,v=self.features.transform(jj);out=contrast_forward(params,c,v,self.data['study'][jj],z['log_positive_weight'],backend='numpy')
                for k in values:values[k].append(out[k])
            result={k:np.concatenate(v) for k,v in values.items()}
        elif family=='conditional_avi_stacker':
            from conditional_stacking import predict_models
            base=self.predict(registry[record['config']['base_model']],ix,registry);expert=base['conditional_scores'][:,1] if 'conditional_scores' in base else base['overall']
            with (OUT/record['artifact']).open('rb') as f:models=pickle.load(f)
            score=predict_models(models,expert,self.avi_features()[ix],sid,record['config']['include_DNA_expert']);result=dict(base);inactive=base['conditional_scores'][:,0] if 'conditional_scores' in base else np.zeros(len(ix));result['conditional_scores']=np.c_[inactive,score]
        else:raise ValueError(family)
        for name,p in result.items():assert len(p)==len(ix) and np.isfinite(p).all(),(record['name'],name)
        self.cache[key]=result;return result

def sigmoid(z):return 1/(1+np.exp(-np.clip(np.asarray(z,dtype='float64'),-30,30)))
def gelu(x):return .5*x*(1+np.tanh(.7978845608028654*(x+.044715*x**3)))
def sha_indices(ix):return hashlib.sha256(np.ascontiguousarray(ix).tobytes()).hexdigest()
