"""Freeze both arms, confidence and all reserved predictions before target join."""
import sys,itertools,time,joblib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
from coverage_validation_common import *
from calibrate_coverage_confidence import uncertainty
from audit_and_lock_library_measurements import write_gzip_csv
SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'

def main():
    assert_previous_sources();lock_path=FINAL/'before_reserved_outcomes_lock.json';assert not lock_path.exists();assert not (FINAL/'reserved_outcome_evaluation.json').exists()
    assert json.loads((FINAL/'fixed_fits_independent_audit.json').read_text())['passed']
    for arm in ['sequence','AVI']:
        audit=json.loads((FINAL/(arm+'_confidence_independent_audit.json')).read_text());assert audit['passed'] and audit['confidence_package_sha256']==sha(FINAL/(arm+'_confidence.joblib'))
    metrics_audit=json.loads((FINAL/'inference_and_metric_pretest_audit.json').read_text());assert metrics_audit['passed']
    for name,h in metrics_audit['code_sha256'].items():assert sha(Path(__file__).parent/name)==h
    receipt=json.loads((DEST/'reserved_input_receipt.json').read_text());assert receipt['passed']
    for name,key in [('reserved_reader_features.npy','features_sha256'),('reserved_physical_metadata.npy','physical_metadata_sha256'),('reserved_observation_metadata.csv.gz','observation_metadata_sha256'),('reserved_physical_sequence_bank.csv.gz','sequence_bank_sha256')]:assert sha(DEST/name)==receipt[key]
    policy=json.loads((FINAL/'confidence_policy_lock.json').read_text());rows=list(csv_rows(DEST/'reserved_observation_metadata.csv.gz'));physical=list(csv_rows(DEST/'reserved_physical_sequence_bank.csv.gz'));banks=np.array([int(r['bank_index']) for r in rows]);bank=np.load(DEST/'reserved_reader_features.npy');pp=np.load(DEST/'reserved_physical_metadata.npy');x=np.c_[bank.reshape(len(bank),48),pp[:,0]/230.,pp[:,1]/np.maximum(pp[:,0]-1,1),np.eye(16)[pp[:,2]*4+pp[:,3]]].astype('float32')
    avi_receipt=json.loads((DEST/'coverage_AVI_input_receipt.json').read_text());assert avi_receipt['passed'];assert sha(DEST/'coverage_AVI_comparison_inputs.csv.gz')==avi_receipt['comparison_inputs_sha256'];amap={r['variant_GRCh38']:r for r in csv_rows(DEST/'coverage_AVI_comparison_inputs.csv.gz')};avi=[]
    for r in rows:
        a=amap[r['variant_GRCh38']];avi.append([float(a['raw']) if a['raw'] else np.nan,float(a['phred']) if a['phred'] else np.nan,int(a['available'])])
    avi=np.asarray(avi,dtype='float32');summaries={}
    for arm in ['sequence','AVI']:
        package=joblib.load(FINAL/(arm+'_confidence.joblib'));assert package['policy_lock_sha256']==sha(FINAL/'confidence_policy_lock.json');meta,models,_=load_ensemble('none',arm)
        cell=np.array([meta['cell_names'].index(r['cell_group']) for r in rows],dtype='int64');known=np.array([r['exact_line_known']=='True' for r in rows]);line=np.array([meta['cell_line_names'].index(r['cell_line']) if r['cell_line'] in meta['cell_line_names'] else meta['cell_line_names'].index('K562') for r in rows],dtype='int64')
        pred=predict_inputs(models,meta,bank[banks],x[banks],cell,line,known,avi);features,b,cm,cs,ls,counts=confidence_inputs(pred,x[banks],cell,line,meta,known);u=uncertainty(package,features,cell,b,cs,ls,counts,policy);assert not np.any(u['pretest_policy_accepted'][~known]);np.savez_compressed(FINAL/(arm+'_reserved_predictions.npz'),row_index=np.arange(len(rows)),**pred,baseline=b.astype('float32'),cell_baseline=cm.astype('float32'),**u)
        summaries[arm]={'rows':len(rows),'pretest_accepted':int(u['pretest_policy_accepted'].sum()),'unknown_lines_accepted':int(u['pretest_policy_accepted'][~known].sum()),'prediction_file_sha256':sha(FINAL/(arm+'_reserved_predictions.npz'))}
        print('RESERVED PREDICTIONS FROZEN',arm,summaries[arm],flush=True)
    cross=[];local=[];same=collections.defaultdict(list)
    for i,r in enumerate(rows):same[(r['source'],r['library'],int(r['bank_index']))].append(i)
    for ids in same.values():
        for a,b in itertools.combinations(ids,2):
            if rows[a]['cell_line']!=rows[b]['cell_line']:cross.append((a,b))
    def pos(i):return int(physical[int(rows[i]['bank_index'])]['position_hg38'])
    order=sorted(range(len(rows)),key=lambda i:(rows[i]['source'],rows[i]['library'],rows[i]['cell_line'],rows[i]['component'],pos(i),rows[i]['variant_GRCh38'],i))
    for a,b in zip(order,order[1:]):
        ra,rb=rows[a],rows[b]
        if all(ra[k]==rb[k] for k in ['source','library','cell_line','component']) and 0<=pos(b)-pos(a)<=25 and ra['variant_GRCh38']!=rb['variant_GRCh38'] and ra['bank_index']!=rb['bank_index']:local.append((a,b))
    np.savez_compressed(FINAL/'reserved_measured_pair_identities.npz',local=np.array(local,dtype='int32').reshape(-1,2),cross=np.array(cross,dtype='int32').reshape(-1,2))
    write_json(FINAL/'reserved_prediction_receipt.json',{'arms':summaries,'pair_counts':{'local':len(local),'cross':len(cross)},'new_numeric_targets_joined':False,'authorized_project_answers':0})
    # Include prior adverse result locks in the immutability check as well.
    previous=ROOT/'outputs/sleep_frailty_ai_v9_validation/before_fresh_outcomes_lock.json';old_lock=json.loads(previous.read_text())
    for path,h in old_lock['files'].items():assert sha(path)==h,path
    files={}
    for path in FINAL.rglob('*'):
        if path.is_file():files[str(path)]=sha(path)
    for path in Path(__file__).parent.glob('*.py'):files[str(path)]=sha(path)
    for name in ['coverage_training_arrays.npz','coverage_measured_pairs.npz','coverage_development_assessment.json','coverage_independent_data_audit.json','library_measurement_role_lock.json','coverage_AVI_input_receipt.json','coverage_AVI_comparison_inputs.csv.gz','reserved_input_receipt.json','reserved_observation_metadata.csv.gz','reserved_physical_sequence_bank.csv.gz','reserved_reader_features.npy','reserved_physical_metadata.npy']:
        files[str(DEST/name)]=sha(DEST/name)
    for path in [SCREEN/'Tardaguila2026_TableS3.tsv',SCREEN/'candidate_strict_no_published_pretraining_neighborhood.csv',SCREEN/'Siraj_unused_libraries/pretraining_guarded_source_roles.csv',previous]:files[str(path)]=sha(path)
    retrieval=json.loads((SCREEN/'Siraj_unused_libraries/targeted_retrieval_receipt.json').read_text())
    for item in retrieval['files']:files[str(ROOT/item['path'])]=sha(ROOT/item['path'])
    write_json(lock_path,{'files':files,'all_reserved_predictions_and_confidence_fixed':True,'reserved_outcome_join_or_performance_calculated':False,'prior_V9_locked_files_verified_unchanged':True,
                        'future_policy':'These cases cease to be untouched after evaluation. Do not tune thresholds, inputs, epochs or model selection against them and call them independent validation again.'})
    print('BEFORE RESERVED OUTCOMES LOCK COMPLETE',sha(lock_path),flush=True)

if __name__=='__main__':main()
