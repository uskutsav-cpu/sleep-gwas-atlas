"""Read selected models, then evaluate the reused benchmark and external assay.

Does not refit, tune, change membership or select by these results.
"""
import csv,json,hashlib,sys
import numpy as np
from common import *
from inference import infer_record

data,seq,avi,manifest,tr,va,te=load();y=data['label'];all_records=json.loads((OUT/'search_results.json').read_text());records={r['name']:r for r in all_records};selection=json.loads((OUT/'model_selection.json').read_text());selected=selection['selected_model']
student_config=json.loads((OUT/'student_config.json').read_text()) if (OUT/'student_config.json').exists() else {}
development_notes=['Round1 benchmark and separate-study performance were already observed before deciding to run the exploratory auxiliary extension. These subsequent evaluations are reused development comparisons, not fresh independent validation.'] if student_config.get('post_round1_extension') else []
assert selected==max(all_records,key=lambda r:r['validation']['average_precision'])['name']
all_features={'avi':avi,'kmer':seq['kmer'],'kmer_avi':np.c_[seq['kmer'],avi],'hyena':seq['hyena_original'],'hyena_avi':np.c_[seq['hyena_original'],avi]}
pre={'avi_training_median_raw':float(np.nanmedian(avi[tr,0])),'avi_training_median_phred':float(np.nanmedian(avi[tr,1])),'missing_avi_training_rows':int(np.isnan(avi[tr,0]).sum()),'selection_sha256':sha(OUT/'model_selection.json'),'evaluation_code_sha256':sha(Path(__file__))}
(OUT/'preprocessing_receipt.json').write_text(json.dumps(pre,indent=2)+'\n')
saved_validation=np.load(OUT/'validation_predictions.npz');replay={};alphabet=np.array(list('ACGT'))
def pairs_from_indices(ix):return [{'variant':str(i),'reference_sequence':''.join(alphabet[data['ref'][i]-7]),'alternate_sequence':''.join(alphabet[data['alt'][i]-7])} for i in ix]
# Numerical prediction replay on held validation data before reporting benchmark.
vfeatures={k:v[va] for k,v in all_features.items()};vpairs=pairs_from_indices(va)
for name,record in records.items():
    actual=infer_record(record,vfeatures,records,vpairs);expected=saved_validation[name];error=float(np.max(np.abs(actual-expected)));ap=metrics(y[va],actual)['average_precision'];difference=abs(ap-record['validation']['average_precision'])
    assert error<3e-5 and difference<1e-4,(name,error,difference)
    replay[name]={'maximum_reload_prediction_error':error,'validation_AP_difference':difference}
(OUT/'checkpoint_replay.json').write_text(json.dumps(replay,indent=2)+'\n')

# Keep every fixed model visible; overall recommendation remains validation-only.
benchmark_features={k:v[te] for k,v in all_features.items()};benchmark_pairs=pairs_from_indices(te);bp={};br={}
for name,record in records.items():
    bp[name]=infer_record(record,benchmark_features,records,benchmark_pairs);br[name]=metrics(y[te],bp[name])
np.savez_compressed(OUT/'reused_benchmark_predictions.npz',row_index=te,label=y[te],**bp)
(OUT/'reused_benchmark_results.json').write_text(json.dumps({'scope':'Previously observed v1 chromosomes 19-22; not a new independent test','selected_by_validation':selected,'models':br},indent=2)+'\n')
print('REUSED BENCHMARK selected',selected,json.dumps(br[selected]),flush=True)

ef=WORK/'external_features_local.npz';ep=WORK/'external_pairs_local.csv';ea=OUT/'external_avi_retrieval_audit.json'
if not all(p.exists() for p in [ef,ep,ea]):
    print('External inputs not complete: benchmark saved, external evaluation pending',flush=True);sys.exit(0)
receipt=json.loads(ea.read_text());assert not receipt['failed_groups'] and not receipt['partial_run'];assert receipt['output_sha256']==sha(OUT/'external_avi_exact_allele_scores.csv')
external=np.load(ef);erows=list(csv.DictReader(ep.open()));assert len(erows)==len(external['label'])
av={r['variant_GRCh38']:r for r in csv.DictReader((OUT/'external_avi_exact_allele_scores.csv').open())};values=[]
for r in erows:
    key=f"{r['chromosome']}:{r['position_hg38']}:{r['genome_reference_allele']}:{r['tested_alternate_allele']}";a=av[key]
    values.append([float(a['avi_raw']),float(a['avi_phred']),1] if a['status']=='exact_allele_match' else [np.nan,np.nan,0])
aa=np.array(values,dtype='float32');features={'avi':aa,'kmer':external['kmer'],'kmer_avi':np.c_[external['kmer'],aa],'hyena':external['hyena_original'],'hyena_avi':np.c_[external['hyena_original'],aa]};ey=external['label'];predictions={};res={}
for name,record in records.items():
    predictions[name]=infer_record(record,features,records,erows);res[name]=metrics(ey,predictions[name])
# No model reselection here, including if another model scores better externally.
result={'scope':'Separate neuronal study, evaluated without fitting its labels; region exclusion and native-context verification applied','selected_by_Salomon_validation':selected,'models':res,'external_rows':len(ey),'external_positive_calls':int(ey.sum()),'limitations':json.loads((OUT/'external_mapping_audit.json').read_text())['limitations']+['AVI training/benchmark overlap is unresolved, so full foundation-model independence is not established','Prior broader-project analyses inspected this study; this is external to current fitting, not a completely unseen publication','This is not a sleep/frailty prediction cohort or genetic-score validation']+development_notes,'selection_sha256':sha(OUT/'model_selection.json')}
(OUT/'external_assay_results.json').write_text(json.dumps(result,indent=2)+'\n')
# Local-only predictions because source article forbids sharing adapted material.
np.savez_compressed(WORK/'external_predictions_local.npz',label=ey,**predictions)
print('EXTERNAL selected',selected,json.dumps(res[selected]),flush=True)
