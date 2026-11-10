"""Build and smoke-test a local inference/metric reconstruction bundle.

This is a private local artifact. It does not publish, upload or claim a fresh
full-training reconstruction or compatibility with another operating system.
"""
import sys,shutil,subprocess,zipfile,importlib.metadata as md
from common import *
DEST=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ai_v3_delivery')
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def main():
    registry=json.loads((OUT/'final_model_registry.json').read_text());report=json.loads((OUT/'final_report_receipt.json').read_text());assert report['registry_sha256']==sha(OUT/'final_model_registry.json');assert report['text_report_sha256']==sha(OUT/'MODEL_REPORT.txt');assert (OUT/'final_figure_receipt.json').is_file();assert (OUT/'frozen_metric_reproduction_check.json').is_file();DEST.mkdir(exist_ok=True);stage=DEST/'sleep_frailty_ai_v3';assert not stage.exists(),'Do not replace a previous bundle silently';stage.mkdir();copied={}
    def copy(path,relative=None):
        path=Path(path);rel=Path(relative) if relative else path.relative_to(ROOT);target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target);assert digest(path)==digest(target);copied[str(rel)]=digest(target)
    for path in WORK.glob('*.py'):copy(path)
    copy(WORK/'requirements_reproduction.txt')
    for name in ['model.py','model_mlx.py']:copy(OLD/name)
    for folder in [OLD/'runtime',OLD/'hyenadna',WORK/'motifs',WORK/'reference_chromatin']:
        for path in folder.rglob('*'):
            if path.is_file() and not path.name.startswith('._') and '__pycache__' not in path.parts:copy(path)
    for name in ['hg38ToHg19.over.chain.gz','hg19ToHg38.over.chain.gz']:copy(ROOT/'work'/name)
    for stem in ['salomon','siraj']:copy(WORK/(stem+'_manifest_local.csv'))
    site=ROOT/'work/stage5_env/lib/python3.13/site-packages'
    for folder in list(site.glob('pyliftover'))+list(site.glob('pyliftover-*.dist-info')):
        for path in folder.rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts:copy(path)
    skip={'model_card_draft.json','progress_snapshot.json','run_handles.json','completion_requirements.json','local_bundle_receipt.json'}
    for path in OUT.iterdir():
        if path.is_file() and path.suffix in ['.json','.csv','.txt','.png','.svg'] and path.name not in skip and not path.name.startswith('._'):copy(path)
    copy(OUT/'reserved_test_predictions.npz')
    for record in registry['records'].values():
        if record['artifact']:
            path=OUT/record['artifact'];assert digest(path)==record['checkpoint_sha256'];copy(path)
    copy(ROOT/'outputs/ALPHAGENOME_OUTPUT_TERMS_NOTICE.md')
    for name in ['avi_exact_allele_scores.csv','avi_retrieval_audit.json']:copy(ROOT/'outputs/sleep_frailty_ai_v2'/name)
    copy(ROOT/'work/sleep_frailty_ml_v2/independent_metrics.R')
    text=['LOCAL MOLECULAR AI INFERENCE AND METRIC RECONSTRUCTION','', 'Read outputs/sleep_frailty_ai_v3/MODEL_REPORT.txt first. This is a research model of laboratory DNA readouts, not a patient sleep/frailty or genetic-score predictor. This local bundle does not publish or upload anything. Source terms and attributions remain in the receipts; this is not a declaration of unrestricted redistribution.','', 'Tested runtime: existing macOS ARM Python3.13.13 with package versions in runtime_metadata.json. The task-local MLX/Autograd and pyliftover files are included. NumPy, SciPy and scikit-learn are supplied by the tested Python environment. A fresh install/another platform is not validated.','', 'From this folder, score a paired-DNA JSON list:',f'python3 work/sleep_frailty_ml_v3/predict.py --input outputs/sleep_frailty_ai_v3/raw_inference_validation_example.json --output example_scores.json --model {registry["role_model_names"]["within_active_variant_ranker"]}','', 'Reconstruct metrics from immutable saved predictions (no new inference or validation):','python3 work/sleep_frailty_ml_v3/reproduce_frozen_metrics.py','', 'The bundle includes trained final checkpoints, frozen predictions, reference/motif/encoder inputs, exact AlphaGenome lookup rows, raw-input examples, evaluation/selection/source audits and all training source code. Large prepared training-feature matrices and the complete source raw datasets are retained in the original workspace/research SSD and identified by hashes in the receipts. This bundle alone is NOT a self-contained fresh retraining package. Use the original worktree and cached inputs for the existing complete training reconstruction. Do not refit/reselect the published V3 artifact after its reserved evaluation.','', 'Input coordinates must be caller-verified GRCh38 substitution coordinates. DNA strings must differ at exactly one base and have200or270bases; assay_context is Salomon_NGN2 or Siraj_OL42_SKNSH. Scores across different contexts are not comparable cell-type relevance or disease-risk measures.'];(stage/'REPRODUCTION.txt').write_text('\n'.join(text)+'\n')
    model=registry['role_model_names']['within_active_variant_ranker'];example='outputs/sleep_frailty_ai_v3/raw_inference_validation_example.json';expected=OUT/'final_cli_example_scores.json';subprocess.run([sys.executable,str(WORK/'predict.py'),'--input',str(ROOT/example),'--output',str(expected),'--model',model],cwd=ROOT,check=True);copy(expected)
    actual=stage/'relocation_example_scores.json';subprocess.run([sys.executable,str(stage/'work/sleep_frailty_ml_v3/predict.py'),'--input',str(stage/example),'--output',str(actual),'--model',model],cwd=stage,check=True);a=json.loads(actual.read_text());e=json.loads(expected.read_text());assert a['input_sha256']==e['input_sha256'] and len(a['scores'])==len(e['scores']);errors=[]
    for x,y in zip(a['scores'],e['scores']):
        assert x['id']==y['id'] and x['variant_GRCh38']==y['variant_GRCh38']
        for k,v in x.items():
            if isinstance(v,(int,float)):errors.append(abs(v-y[k]))
    assert errors and max(errors)<1e-4;subprocess.run([sys.executable,str(stage/'work/sleep_frailty_ml_v3/reproduce_frozen_metrics.py')],cwd=stage,check=True);check={'relocation_max_numeric_score_error':max(errors),'rows':len(a['scores']),'runtime':'Same tested interpreter/environment on macOS ARM; scorer assets and weights relocated to the local bundle','not_proved':'Fresh dependency installation, another operating system or full retraining without the original cached training inputs'};(stage/'relocation_check.json').write_text(json.dumps(check,indent=2)+'\n')
    files=[p for p in stage.rglob('*') if p.is_file() and not p.name.startswith('._') and '__pycache__' not in p.parts];manifest={str(p.relative_to(stage)):digest(p) for p in files};(stage/'FILE_SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n');archive=DEST/'sleep_frailty_ai_v3_local_bundle.zip';assert not archive.exists()
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path in files+[stage/'FILE_SHA256.json']:z.write(path,Path(stage.name)/path.relative_to(stage))
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None;assert len(z.infolist())==len(files)+1
    output={'archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_sha256':digest(archive),'staged_directory':str(stage),'file_count':len(files)+1,'registry_sha256':digest(OUT/'final_model_registry.json'),'report_receipt_sha256':digest(OUT/'final_report_receipt.json'),'relocation_check':check,'zip_CRC_check_passed':True,'scope':'Local frozen-model inference and metric reconstruction, with training code/provenance. Large training matrices/raw inputs remain in the original worktree/SSD; fresh retraining elsewhere not verified.'};(OUT/'local_bundle_receipt.json').write_text(json.dumps(output,indent=2)+'\n');print('LOCAL VERIFIED BUNDLE',archive,'bytes',archive.stat().st_size,flush=True)
if __name__=='__main__':main()
