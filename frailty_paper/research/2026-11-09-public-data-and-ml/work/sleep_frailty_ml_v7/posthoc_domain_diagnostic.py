"""Diagnose failed transfer after tests; never alter locks/models/thresholds."""
import gzip,collections
from common import *
from sequence_features import Encoder
def main():
    assert (OUT/'reserved_study_results.json').exists() and (OUT/'additional_neural_results.json').exists();lock=json.loads((OUT/'final_model_and_prediction_lock.json').read_text());record=lock['sequence_only']['records'][0]['record']
    with np.load(record['transform']) as z:mean=z['mean'];sd=z['sd']
    report={'posthoc_after_outside_performance':True,'no_models_thresholds_or_original_results_changed':True,'datasets':{}}
    for label,path in [('primary_native150',DATA/'reserved_feature_bank.npy'),('additional_native199_context151',DATA/'additional_neural_feature_bank.npy')]:
        x=np.load(path,mmap_mode='r');z=abs((np.asarray(x)-mean)/sd);maximum=z.max(1);top=np.argmax(z,axis=1);report['datasets'][label]={'n':len(x),'fraction_above_domain_limit12':float(np.mean(maximum>12)),'maximum_z_quantiles':np.quantile(maximum,[0,.5,.9,1]).tolist(),'most_frequent_maximum_columns':dict(collections.Counter(map(int,top)).most_common(8))}
    old=np.load(DATA/'feature_bank.npy',mmap_mode='r');ids=np.linspace(0,len(old)-1,1000,dtype=int);z=abs((np.asarray(old[ids])-mean)/sd);report['datasets']['evenly_spaced_development_inputs']={'n':len(ids),'fraction_above_domain_limit12':float(np.mean(z.max(1)>12)),'maximum_z_quantiles':np.quantile(z.max(1),[0,.5,.9,1]).tolist()}
    examples=[]
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:
        for r in csv.DictReader(f):
            p=int(r['changed_position']);a,b=r['reference'],r['alternate'];start=p-75
            if start>=0 and start+151<=len(a):examples.append((a[start:start+151],b[start:start+151],75))
            if len(examples)==6:break
    encoder=Encoder();full=encoder.features([a for a,b,p in examples],[b for a,b,p in examples],[p for a,b,p in examples]);short=encoder.features([a[:-1] for a,b,p in examples],[b[:-1] for a,b,p in examples],[p for a,b,p in examples]);fullz=abs((full-mean)/sd);shortz=abs((short-mean)/sd);report['terminal_base_crop_control']={'full_length':151,'cropped_length':150,'SNP_and_variant51bp_sequence_identical':True,'maximum_change_outside_native_context_columns':float(np.max(abs(full[:,128:]-short[:,128:]))),'native_context_feature_max_changes':np.max(abs(full[:,:128]-short[:,:128]),axis=1).tolist(),'full_maximum_z':fullz.max(1).tolist(),'cropped_maximum_z':shortz.max(1).tolist(),'cropped_maximum_columns':np.argmax(shortz,axis=1).tolist(),'scope':'Numerical input-length diagnostic on development DNA, not a biological experiment or a new model-performance test.'};write_json(OUT/'posthoc_context_length_diagnostic.json',report)
    joined=list(csv.DictReader((OUT/'reserved_all_observations.csv').open()));by={(int(r['design_index']),r['cell']):r for r in joined};ids=sorted({i for i,c in by if c=='HNP'}&{i for i,c in by if c=='HEK293'});actual=np.array([float(by[(i,'HNP')]['measured_absolute_log2_effect'])-float(by[(i,'HEK293')]['measured_absolute_log2_effect']) for i in ids]);zero=float(np.mean(abs(actual)));result={'posthoc_descriptive_null_comparator':True,'matched_variants':len(ids),'zero_between_cell_difference_MAE':zero,'arms':{}}
    for arm in ['sequence_only','AVI_added']:
        p=np.array([float(by[(i,'HNP')][arm+'_predicted_absolute_log2_effect'])-float(by[(i,'HEK293')][arm+'_predicted_absolute_log2_effect']) for i in ids]);mae=float(np.mean(abs(actual-p)));result['arms'][arm]={'MAE':mae,'MAE_skill_vs_zero_cell_difference':1-mae/zero,'fraction_of_predicted_differences_positive':float(np.mean(p>0))}
    result['scope']='Added after unblinding to avoid interpreting improvement over a shifted training-cell median as learned variant-specific cell interactions. No refitting, retention-rule change or independent confirmatory claim.';write_json(OUT/'posthoc_zero_cell_difference_comparator.json',result);print('POSTHOC DOMAIN AND CELL-DIFFERENCE DIAGNOSTICS WRITTEN',flush=True)
if __name__=='__main__':main()
