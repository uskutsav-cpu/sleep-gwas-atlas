"""Build and locally relocate-check a compact inference delivery bundle."""
import shutil,subprocess,sys,zipfile,importlib.metadata,os
from common import *
sys.path.insert(0,str(OLD/'runtime'))
def main():
    assert (OUT/'reserved_study_results.json').exists();delivery=DATA.parent/'sleep_frailty_ai_v6_delivery';delivery.mkdir(exist_ok=True);stage=delivery/'sleep_frailty_ai_v6';stage.mkdir(exist_ok=True);bw=stage/'work/sleep_frailty_ml_v6';bo=stage/'model_data/outputs';bm=stage/'model_data/models';bw.mkdir(parents=True,exist_ok=True);bo.mkdir(parents=True,exist_ok=True);bm.mkdir(parents=True,exist_ok=True)
    names=['common.py','data.py','network.py','predict_checkpoints.py','prepare_features.py','encode_inputs.py','calibrate_uncertainty.py','score.py']
    for name in names:shutil.copyfile(WORK/name,bw/name)
    text=(bw/'common.py').read_text();text=text.replace("DATA=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v6')","DATA=ROOT/'model_data'");assert "DATA=ROOT/'model_data'" in text;(bw/'common.py').write_text(text)
    text=(bw/'predict_checkpoints.py').read_text();text=text.replace("path=Path(record['checkpoint']);assert sha(path)==record['checkpoint_sha256']","path=Path(record['checkpoint']);path=path if path.is_absolute() else ROOT/path;assert sha(path)==record['checkpoint_sha256']");assert 'else ROOT/path' in text;(bw/'predict_checkpoints.py').write_text(text)
    registry=json.loads((OUT/'final_fit_registry.json').read_text());policy=json.loads((OUT/'uncertainty_policy.json').read_text())
    for r in registry['records']:
        source=Path(r['checkpoint']);target=bm/source.name;shutil.copyfile(source,target);assert sha(target)==r['checkpoint_sha256'];r['checkpoint']='model_data/models/'+source.name
    write_json(bo/'final_fit_registry.json',registry);policy['primary_registry_sha256']=sha(bo/'final_fit_registry.json');write_json(bo/'uncertainty_policy.json',policy)
    for name in ['fixed_features_receipt.json','selected_motifs.json','PLAIN_LANGUAGE_RESULTS.txt','model_card.json','readiness_assessment.json','reserved_study_results.json','reserved_block_bootstrap.json','reserved_additional_checks.json','model_search_plan.json','development_selection.json','partition_audit.json','blood_extension_decision.json','immutable_reuse_audit.json','K562_source_allele_audit.json','final_checkpoint_audit.json','scorer_smoke_audit.json','uncertainty_protocol_details.json','kosicki_continuous_results.json','kosicki_continuous_recovery_audit.json','reserved_source_missingness_diagnostic.json','additional_study_blind_design_audit.json','additional_reference_audit.json','final_model_lock.json','additional_study_validation_lock.json','continuous_recovery_evaluation_lock.json','additional_analysis_code_lock.json','additional_protocol_extension_receipt.json','additional_prior_exposure_search.json','neuronal_posthoc_measurement_reliability.json']:
        shutil.copyfile(OUT/name,bo/name)
    old=stage/'work/sleep_frailty_ml';old.mkdir(parents=True,exist_ok=True)
    for name in ['model.py','model_mlx.py']:shutil.copyfile(OLD/name,old/name)
    (old/'hyenadna').mkdir(exist_ok=True);shutil.copyfile(OLD/'hyenadna/model.safetensors',old/'hyenadna/model.safetensors')
    for name in ['config.json','README.md','LICENSE','LICENSE.txt']:
        if (OLD/'hyenadna'/name).exists():shutil.copyfile(OLD/'hyenadna'/name,old/'hyenadna'/name)
    runtime=old/'runtime';runtime.mkdir(exist_ok=True)
    for name in ['mlx','mlx-0.22.1.dist-info','autograd','autograd-1.8.0.dist-info']:
        source=OLD/'runtime'/name
        if source.exists():shutil.copytree(source,runtime/name,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','._*'))
    motif=stage/'work/sleep_frailty_ml_v3';(motif/'motifs').mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/'work/sleep_frailty_ml_v3/motif_model.py',motif/'motif_model.py')
    for name in ['H13CORE_annotation.jsonl','H13CORE_pwms.txt']:shutil.copyfile(ROOT/'work/sleep_frailty_ml_v3/motifs'/name,motif/'motifs'/name)
    provenance=stage/'training_source';provenance.mkdir(exist_ok=True)
    for p in WORK.iterdir():
        if p.is_file() and not p.name.startswith('._') and p.suffix in ['.py','.R','.txt']:shutil.copyfile(p,provenance/p.name)
    examples=stage/'examples';examples.mkdir(exist_ok=True);shutil.copyfile(DATA/'scorer_smoke_inputs.csv',examples/'pairs.csv');shutil.copyfile(OUT/'scorer_smoke_outputs.csv',examples/'expected_scores.csv')
    requirements={name:importlib.metadata.version(name) for name in ['numpy','scipy','scikit-learn','mlx','autograd']};(stage/'requirements.txt').write_text('\n'.join(name+'=='+version for name,version in requirements.items())+'\n')
    readme=['SLEEP–FRAILTY RESEARCH DNA-EFFECT MODEL — LOCAL INFERENCE BUNDLE','',
        'This experimental bundle ranks paired DNA changes and can leave weak/unsupported inputs unranked. It failed the neuronal effect-ranking check and its confidence policy failed the external liver-study check. Ranked statuses and probabilities do not have validated cross-study confidence; broad biological reliability is not established. Read both external study results before using it. It does not predict a person’s insomnia/frailty risk or establish causal genes.',
        '', 'RUN', 'Use Python 3.13 on an Apple Silicon Mac with the packages in requirements.txt. MLX/autograd runtime files are included from the tested local installation; NumPy/SciPy/scikit-learn are external dependencies.',
        'python3 work/sleep_frailty_ml_v6/score.py --input examples/pairs.csv --output my_scores.csv',
        'Run that command from the bundle directory, or use absolute paths. Input columns:id,reference,alternate,cell_group. A blank ranked_score means the model abstained; do not replace it with zero.',
        '', 'WHAT WAS TESTED', 'The staged bundle was executed separately from the original working tree on the same Mac/interpreter and its outputs were compared with the original scorer. This is a local relocation check, not a fresh-platform installation test.',
        'Checkpoint bytes and uncertainty thresholds are unchanged. Only portable file paths and the corresponding registry/policy metadata hash were rebased. Original frozen scientific reports and model hashes remain in model_data/outputs.',
        '', 'WHAT IS INCLUDED','Final model checkpoints; fixed DNA encoder; motif resources; inference source; example inputs/expected outputs; scientific reports; data and model receipts; full V6 source snapshot for provenance.',
        'Large training caches and raw assay archives remain in the original external-SSD workspace. The provenance source snapshot is not a one-command, fresh-platform full-retraining package.',
        '', 'PUBLIC ATTRIBUTION','Kosicki et al.2025,doi10.1038/s41467-025-60064-1; Zenodo15353561; CC BY4.0.', 'IGVFDS6149YGOL HepG2 measurement study / IGVFDS1589ODOW library; current released IGVF MPRA Quantification v0.0.7; https://api.data.igvf.org/analysis-sets/IGVFDS6149YGOL/?format=json.', 'Keukeleire et al.2025,doi10.1186/s12859-025-06065-9; BCalm-associated public assay; actual current-release sample/denominators documented in the model card.', 'Kircher et al.2019,doi10.1038/s41467-019-11526-w; saturation-MPRA portal; article CC BY4.0; portal codeGPL3.', 'Siraj et al.2026,doi10.1038/s41586-026-10121-6; Zenodo15297965.', 'Salomon et al.2026,PMC13405209; kircherlab/80K-Analysis.', 'Original tiny HyenaDNA checkpoint revisione8c1effa8673814e257e627d2e1eda9ea5a373f6; checkpointSHA2565ce2146c21e9c4baa6bddc4998fd3d029903ae84a563bf80218644082194a12d.', 'HOCOMOCOv13human motif resource, declared subset in selected_motifs.json. Research source/data licences and original limitations continue to apply.',
        '', 'READ FIRST','model_data/outputs/PLAIN_LANGUAGE_RESULTS.txt; model_card.json; readiness_assessment.json.',
    ];(stage/'README.txt').write_text('\n'.join(readme)+'\n')
    result=stage/'examples/relocated_scores.csv';command=[sys.executable,str(bw/'score.py'),'--input',str(examples/'pairs.csv'),'--output',str(result)];run_env=os.environ.copy();run_env['PYTHONPATH']=os.pathsep.join(str(Path(v).resolve()) for v in run_env.get('PYTHONPATH','').split(os.pathsep) if v);subprocess.run(command,check=True,cwd=str(stage),env=run_env);expected=list(csv.DictReader((examples/'expected_scores.csv').open()));got=list(csv.DictReader(result.open()));assert len(expected)==len(got);maximum=0.
    for a,b in zip(expected,got):
        assert a['id']==b['id'] and a['status']==b['status'] and a['prediction_set']==b['prediction_set']
        for name in ['ranked_score','calibrated_laboratory_call_probability']:
            if not a[name]:assert not b[name]
            else:maximum=max(maximum,abs(float(a[name])-float(b[name])))
    assert maximum<1e-6,maximum
    files={str(p.relative_to(stage)):sha(p) for p in stage.rglob('*') if p.is_file() and not any(part.startswith('._') or part=='__pycache__' for part in p.relative_to(stage).parts)};write_json(stage/'FILE_SHA256.json',files);archive=delivery/'sleep_frailty_ai_v6_local_inference.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name in sorted(files):z.write(stage/name,'sleep_frailty_ai_v6/'+name)
        z.write(stage/'FILE_SHA256.json','sleep_frailty_ai_v6/FILE_SHA256.json')
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    receipt={'archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':sha(archive),'stage_directory':str(stage),'file_count':len(files),'all_ZIP_CRCs_passed':True,'staged_model_checkpoint_bytes_match_original':True,'actual_relocated_scorer_invoked':True,'maximum_output_probability_difference':maximum,'local_interpreter':sys.executable,'packages':requirements,'scope':'Inference bundle, locally relocated on the same testedMac/interpreter. Not independent biology, a new held-out test, fresh-platform portability or full retraining reconstruction.'};write_json(OUT/'local_delivery_receipt.json',receipt);print('LOCAL INFERENCE DELIVERY BUILT AND REPLAYED',archive.stat().st_size,maximum,flush=True)
if __name__=='__main__':main()
