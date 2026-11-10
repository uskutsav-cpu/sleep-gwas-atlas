"""Build a locally replayed inference bundle without changing trained weights."""
import shutil,subprocess,zipfile,importlib.metadata,sys
from common import *
sys.path.insert(0,str(OLD/'runtime'))
def main():
    assert json.loads((OUT/'final_checkpoint_and_R_audit.json').read_text())['passed'];delivery=DATA.parent/'sleep_frailty_ai_v7_delivery';stage=delivery/'sleep_frailty_ai_v7';bw=stage/'work/sleep_frailty_ml_v7';bo=stage/'model_data/outputs';bm=stage/'model_data/models'
    for p in [bw,bo,bm]:p.mkdir(parents=True,exist_ok=True)
    for name in ['common.py','sequence_features.py','data.py','network.py','train.py','metrics.py','baselines.py','finalize.py','score.py']:shutil.copyfile(WORK/name,bw/name)
    p=bw/'common.py';text=p.read_text().replace("DATA=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v7')","DATA=ROOT/'model_data'").replace("V6=DATA.parent/'sleep_frailty_ml_v6'","V6=ROOT/'reference_v6'");assert "DATA=ROOT/'model_data'" in text and "V6=ROOT/'reference_v6'" in text;p.write_text(text)
    registry=json.loads((OUT/'inference_registry.json').read_text())
    for member in registry['records']:
        r=member['record']
        for key in ['checkpoint','transform']:
            source=Path(r[key]);shutil.copyfile(source,bm/source.name);assert sha(bm/source.name)==r[key+'_sha256'];r[key]='model_data/models/'+source.name
            if key=='transform':shutil.copyfile(source.with_suffix('.json'),bm/source.with_suffix('.json').name)
    write_json(bo/'inference_registry.json',registry)
    for name in ['model_search_plan.json','model_search_plan_initial.json','development_loss_amendment.json','raw_MAE_gradient_audit.json','development_selection_round1.json','development_selection.json','reserved_study_results.json','reserved_all_observations.csv','final_model_and_prediction_lock.json','feature_receipt.json','reserved_design_and_original_source_audit.json','reserved_input_arithmetic_audit.json','development_data_audit.json','compact_data_and_pair_receipt.json','independent_partition_audit.json','AVI_collection_receipt.json','AVI_permissive_artifact_receipt.json','final_checkpoint_and_R_audit.json','PLAIN_LANGUAGE_RESULTS.txt','model_card.json','readiness_assessment.json','project_molecular_predictions.csv','project_integration_audit.json','figure_receipt.json','additional_neural_validation_protocol.json','additional_neural_prediction_lock.json','additional_neural_results.json','additional_neural_observations.csv','combined_model_recommendation.json','neuroimmune_candidate_source_audit.json','extra_AVI_inference_receipt.json','all_project_variant_input_audit.json','complete_project_model_integration_audit.json','all_project_variant_predictions.csv','functional_candidate_predictions.csv','complete_gene_evidence_with_model_annotations.csv']:
        shutil.copyfile(OUT/name,bo/name)
    for name in ['scientific_failure_assessment.json','prioritized_model_improvement_plan.json','scientific_report_refinement_receipt.json','posthoc_context_length_diagnostic.json','posthoc_zero_cell_difference_comparator.json','Matoba_fresh_source_schema_audit.json','Matoba_official_GEO_schema_audit.json']:
        shutil.copyfile(OUT/name,bo/name)
    for path in OUT.glob('development_*_curve.json'):shutil.copyfile(path,bo/path.name)
    if (OUT/'nearby_credible_variant_comparisons.csv').exists():shutil.copyfile(OUT/'nearby_credible_variant_comparisons.csv',bo/'nearby_credible_variant_comparisons.csv')
    shutil.copytree(OUT/'figures',bo/'figures',dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
    if (DATA/'AVI_inference_lookup.csv.gz').exists():shutil.copyfile(DATA/'AVI_inference_lookup.csv.gz',stage/'model_data/AVI_inference_lookup.csv.gz')
    rv=stage/'reference_v6/outputs';rv.mkdir(parents=True,exist_ok=True)
    for name in ['fixed_features_receipt.json','selected_motifs.json']:shutil.copyfile(V6OUT/name,rv/name)
    old=stage/'work/sleep_frailty_ml';(old/'hyenadna').mkdir(parents=True,exist_ok=True)
    for name in ['model.py','model_mlx.py']:shutil.copyfile(OLD/name,old/name)
    for name in ['model.safetensors','config.json','README.md','LICENSE','LICENSE.txt']:
        if (OLD/'hyenadna'/name).exists():shutil.copyfile(OLD/'hyenadna'/name,old/'hyenadna'/name)
    runtime=old/'runtime';runtime.mkdir(exist_ok=True)
    for name in ['mlx','mlx-0.22.1.dist-info','autograd','autograd-1.8.0.dist-info']:
        if (OLD/'runtime'/name).exists():shutil.copytree(OLD/'runtime'/name,runtime/name,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','._*'))
    motif=stage/'work/sleep_frailty_ml_v3';(motif/'motifs').mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/'work/sleep_frailty_ml_v3/motif_model.py',motif/'motif_model.py');shutil.copyfile(ROOT/'work/sleep_frailty_ml_v3/motifs/H13CORE_pwms.txt',motif/'motifs/H13CORE_pwms.txt')
    source=stage/'training_source';source.mkdir(exist_ok=True)
    for p in WORK.iterdir():
        if p.is_file() and p.suffix in ['.py','.R','.txt']:shutil.copyfile(p,source/p.name)
    examples=stage/'examples';examples.mkdir(exist_ok=True);shutil.copyfile(DATA/'scorer_smoke_inputs.csv',examples/'pairs.csv');shutil.copyfile(DATA/'project_model_inputs.csv',examples/'project_pairs.csv');shutil.copyfile(OUT/'scorer_smoke_outputs.csv',examples/'expected.csv');env=environment();env['PYTHONPATH']=os.pathsep.join(str(Path(v).resolve()) for v in env.get('PYTHONPATH','').split(os.pathsep) if v);relocated=examples/'relocated.csv';subprocess.run([sys.executable,str(bw/'score.py'),'--input',str(examples/'pairs.csv'),'--output',str(relocated)],cwd=stage,env=env,check=True)
    expected=list(csv.DictReader((examples/'expected.csv').open()));got=list(csv.DictReader(relocated.open()));assert len(expected)==len(got);maximum=0.
    for a,b in zip(expected,got):
        for k in ['id','status','unranked_reason','input_arm','AVI_available','input_policy_eligible']:assert a[k]==b[k],(k,a[k],b[k])
        for k in ['predicted_absolute_log2_effect','lower90','upper90','ranked_effect_magnitude','raw_model_projection','raw_interval_lower90','raw_interval_upper90']:
            if a[k]=='':assert b[k]==''
            else:maximum=max(maximum,abs(float(a[k])-float(b[k])))
    assert maximum<1e-5
    requirements={name:importlib.metadata.version(name) for name in ['numpy','scipy','scikit-learn','mlx','autograd','joblib']};(stage/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in requirements.items())+'\n')
    (stage/'README.txt').write_text('SLEEP–FRAILTY V7: EXPERIMENTAL DNA-EFFECT MAGNITUDE MODEL\n\nRead model_data/outputs/PLAIN_LANGUAGE_RESULTS.txt and reserved_study_results.json first. A working scorer does not establish useful biological prediction. Actual outside-study errors, coverage and abstention determine what is supported. This is not a patient sleep/frailty risk model and does not establish causal genes.\n\nRun with Python3.13 on Apple Silicon and the listed dependencies:\npython3 work/sleep_frailty_ml_v7/score.py --input examples/pairs.csv --output scores.csv\n\nInput columns: id, reference, alternate, cell_group; optional cell_line and variant_GRCh38. Equal DNA sequences51–1024bp with one unambiguous SNV are required. Unknown cells or unsupported cases remain unranked. A blank ranked_effect_magnitude is not a zero effect. The90%interval label refers to a calibration target; read actual outside coverage before using it.\n\nIf AVI is retained, a score is used only for an exact cached GRCh38 allele and matching51bp reference context. Unavailable values are median-imputed with a missing flag, never declared zero biological effects. AVI is a composite including other predictors.\n\nThe staged scorer was actually relocated and run on the same Mac/interpreter. This is not a fresh-platform installation test. Large training sources remain in the external research workspace. Full source code, data receipts, model-selection rules and scientific reports are included; this is a local inference package, not a one-command complete-retraining release.\n\nSources: McAfee et al.2023,doi10.1016/j.xgen.2023.100404; Koesterich et al.2023,doi10.3390/ijms24043509; IGVFDS9175DMRR and IGVFDS6149YGOL,IGVF240K; source studies and limitations are enumerated in the audits. Prior-development sources include Kosicki2025,Kircher2019,Salomon and Siraj. TinyHyenaDNAcheckpoint e8c1effa (unchanged); HOCOMOCOv13 motifs; public Permissive Use Downloadable AVI SNV scores. Original source licences and research limitations apply.\n')
    (stage/'README.txt').write_text('''SLEEP–FRAILTY V7: FROZEN ADVERSE MOLECULAR BENCHMARK

The outside tests did not establish useful effect-magnitude prediction or reliable variant-specific cell differences. No cell group passed all checks for released rankings. This is not a patient sleep/frailty risk model or evidence of causal genes. Read model_data/outputs/PLAIN_LANGUAGE_RESULTS.txt, scientific_failure_assessment.json and both outside-study result files first.

Small neural heads were trained on measured outcomes. The pretrained DNA reader stayed fixed. AVI reduced absolute error in both studies, but failed the second study's prespecified nearby-ranking gate and was not retained. AVI combines AlphaGenome and other information; its contribution cannot be assigned uniquely to AlphaGenome.

Local replay on the tested Apple Silicon Mac, Python 3.13 and listed dependencies:
python3 work/sleep_frailty_ml_v7/score.py --input examples/pairs.csv --output scores.csv

Input columns: id, reference, alternate, cell_group; optional cell_line, variant_GRCh38 and condition. Equal unambiguous DNA pairs of 51–1024 bases with one single-base substitution are required. Unknown cells, invalid DNA and unsupported stimulation return no prediction. This release withholds authorized prediction and rank fields for every cell group. A blank answer is not a zero or harmless effect. Raw model projections and intervals, when present, are separately labeled unvalidated diagnostics. Broad intervals and withholding all answers do not establish useful confidence.

The scorer was actually relocated and replayed on the same Mac and interpreter. This is not a fresh-platform installation test or a one-command full-retraining release. Original training sources, source audits and frozen predictions remain in the research workspace; source code, receipts, weights and scientific reports are included here. Prior genetic and NHANES findings were preserved.

Primary sources: McAfee et al. 2023, doi:10.1016/j.xgen.2023.100404; Retallick-Townsley et al. 2026, doi:10.1038/s41467-026-77114-x. Additional source identities and limitations are enumerated in the audits. The unchanged TinyHyenaDNA checkpoint, HOCOMOCO motifs, public Permissive Use Downloadable AVI composite and included third-party dependencies retain their original licences. Future-model candidates and the newly accessible Matoba/GSE289858 source are proposed work, not completed validation or demonstrated improvements.
''')
    files={str(p.relative_to(stage)):sha(p) for p in stage.rglob('*') if p.is_file() and '__pycache__' not in p.parts and not p.name.startswith('._') and p.name!='FILE_SHA256.json'};write_json(stage/'FILE_SHA256.json',files);archive=delivery/'sleep_frailty_ai_v7_local_inference.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name in sorted(files):z.write(stage/name,'sleep_frailty_ai_v7/'+name)
        z.write(stage/'FILE_SHA256.json','sleep_frailty_ai_v7/FILE_SHA256.json')
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None;manifest=json.loads(z.read('sleep_frailty_ai_v7/FILE_SHA256.json'))
        for name,h in manifest.items():assert hashlib.sha256(z.read('sleep_frailty_ai_v7/'+name)).hexdigest()==h,name
    write_json(OUT/'local_delivery_receipt.json',{'archive':str(archive),'archive_sha256':sha(archive),'bytes':archive.stat().st_size,'file_count':len(files),'actual_relocated_scorer_executed':True,'maximum_score_difference':maximum,'all_zip_CRCs_passed':True,'scope':'Local inference replay on the same tested Mac and interpreter, not external biological validation or fresh-platform portability.'});print('V7 DELIVERY REPLAYED',archive.stat().st_size,maximum,flush=True)
if __name__=='__main__':main()
