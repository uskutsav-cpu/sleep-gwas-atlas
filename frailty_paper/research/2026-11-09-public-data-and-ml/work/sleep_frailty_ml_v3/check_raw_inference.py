"""Recompute validation-row features from raw DNA and compare model scores.

Only chromosomes4–5 sampled. No reserved-model prediction or outcome is read.
"""
import json,numpy as np
from common import *
from raw_features import RawDNAInference,PublicReferenceChromatin,describe
from inference_records import StoredFeatureInference
from validation_registry import completed_candidates

def main():
    stored=StoredFeatureInference();chosen=[];records=[];offset=0
    for stem in ['salomon','siraj']:
        with np.load(WORK/(stem+'_sequences_local.npz')) as z:a={k:z[k] for k in ['ref','alt','position','chromosome']}
        va=np.flatnonzero(np.isin(a['chromosome'],[4,5]));indices=[int(va[np.argmin(a['position'][va])]),int(va[np.argmax(a['position'][va])])]
        if indices[0]==indices[1]:indices[1]=int(va[-1])
        assert len(set(indices))==2
        for i in indices:
            r=stored.manifest[offset+i];ref=''.join('ACGT'[x-7] for x in a['ref'][i]);alt=''.join('ACGT'[x-7] for x in a['alt'][i]);variant=r['variant_GRCh38'];_,_,gref,galt=variant.split(':');pos=int(a['position'][i]);strand='+' if (ref[pos],alt[pos])==(gref,galt) else '-'
            records.append({'id':'validation_replay_'+str(offset+i),'ref_sequence':ref,'alt_sequence':alt,'chromosome':r['chromosome'],'position_hg38':int(r['position_hg38']),'variant_GRCh38':variant,'genomic_strand':strand,'assay_context':r['study']});chosen.append(offset+i)
        offset+=len(a['chromosome'])
    fixture=(np.array([10,15,40]),np.array([20,25,50]),np.array([20,25,50]),np.array([1.,2.,3.]))
    assert describe(fixture,20,True)[1]==1 and describe(fixture,25,True)[1]==0;assert abs(describe(fixture,30,True)[4]-25/131)<1e-12
    raw=RawDNAInference(records);ix=np.array(chosen);jj=np.arange(len(records));checks={}
    for kind,tolerance in [('dna_pattern',1e-7),('hyena_original',3e-5),('motif',.0006),('chromatin',1e-6),('avi',1e-7)]:
        expected=stored.matrix(kind,ix);actual=raw.matrix(kind,jj);assert np.array_equal(np.isnan(expected),np.isnan(actual));err=float(np.nanmax(np.abs(actual-expected)));assert err<tolerance,(kind,err);checks[kind]={'max_absolute_error':err,'tolerance':tolerance}
    registry,_=completed_candidates();names=['pattern_boost_none_depth3','motif_chromatin_prior_ensemble_activity_effect','motif_chromatin_full_avi_input_ensemble_activity_effect','contrast_prior_val4-5_ensemble_effect','contrast_full_val4-5_ensemble_effect_avi'];predchecks={}
    for name in names:
        expected=stored.predict(registry[name],ix,registry);actual=raw.predict(registry[name],jj,registry);errors={k:float(np.max(np.abs(actual[k]-expected[k]))) for k in expected};assert max(errors.values())<1e-4,(name,errors);predchecks[name]=errors
    p=OUT/'raw_inference_validation_example.json';p.write_text(json.dumps(records,indent=2)+'\n');(OUT/'raw_inference_check.json').write_text(json.dumps({'rows':chosen,'chromosomes':stored.data['chromosome'][ix].tolist(),'feature_checks':checks,'model_checks':predchecks,'scope':'Raw DNA/coordinate pipeline agrees with stored validation feature/checkpoint scores. Not biological validation. No reserved chromosomes1–3 model scores evaluated.','blind_test_evaluated':False},indent=2)+'\n');print('RAW INFERENCE REPLAY PASSED',len(records),'rows',len(names),'models',flush=True)
if __name__=='__main__':main()
