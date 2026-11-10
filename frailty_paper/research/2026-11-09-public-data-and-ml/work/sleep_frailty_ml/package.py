from pathlib import Path
import zipfile,json,hashlib
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v1'
required=['README.md','test_results.json','model_selection.json','inference_check.json','test_uncertainty.json','model_comparison.csv','project_candidate_predictions.csv','project_candidate_embeddings.npz','models/finetuned_hyenadna.npz','models/frozen_hyenadna.npz','models/cnn.npz','models/kmer_baseline.npz']
assert all((OUT/p).is_file() for p in required),[p for p in required if not (OUT/p).is_file()]
start='''# Start here: Sleep–Frailty AI v1

This package contains trained molecular research models, an inference tool, learned variant representations, full benchmark results, and reproducible training code. Read `outputs/sleep_frailty_ai_v1/README.md` for actual results and limitations. It predicts a laboratory neuron-assay call, not a person's insomnia or frailty.

For portable CPU predictions (Python 3.13 tested):

    python3 -m venv .venv
    .venv/bin/pip install numpy==2.4.3 autograd==1.8.0
    .venv/bin/python work/sleep_frailty_ml/predict.py --backend numpy --pairs work/sleep_frailty_ml/project_sequence_pairs.csv --checkpoint outputs/sleep_frailty_ai_v1/models/finetuned_hyenadna.npz --output new_pair_predictions.csv

The same predictor also accepts the frozen checkpoint, the CNN checkpoint and the k-mer baseline checkpoint. Use `model_selection.json` for the validation-based recommendation; the command above illustrates the fine-tuned DNA model. Scores are uncalibrated research outputs, not health-risk probabilities. The 33 example DNA pairs are public reference contexts, not participant genotypes.

To reproduce training on an Apple-silicon Mac, additionally install `mlx==0.22.1`, then run these commands from this extracted root:

    .venv/bin/python work/sleep_frailty_ml/download_inputs.py
    .venv/bin/python work/sleep_frailty_ml/prepare.py
    .venv/bin/python work/sleep_frailty_ml/check_port.py
    .venv/bin/python work/sleep_frailty_ml/train.py
    .venv/bin/python work/sleep_frailty_ml/train_cnn.py
    .venv/bin/python work/sleep_frailty_ml/merge_cnn_results.py
    .venv/bin/python work/sleep_frailty_ml/select_model.py
    .venv/bin/python work/sleep_frailty_ml/evaluate_uncertainty.py
    .venv/bin/python work/sleep_frailty_ml/check_inference.py

Training uses a local GPU and can take substantial time. These scripts recreate fixed experiments; changing models after inspecting test outcomes invalidates a fresh holdout claim. Source inputs are publicly downloadable from pinned revisions and checked by SHA-256. The full training assay table and runtime packages are deliberately omitted from this distribution. Human-reference genome exposure during pretraining and single-assay selection remain limitations.

`EXPERIMENT_PLAN.md` and `CNN_COMPARATOR_PLAN.md` preserve the initial plan and validation-driven debugging amendments. `interrupted_feature_extraction.log` and the unnormalized CNN attempt are retained to make the development history visible. The test results and artifacts include negative findings.

Attribution: LongSafari/HazyResearch HyenaDNA model contributors; Salomon and colleagues' neuron reporter study (DOI 10.64898/2026.07.16.738760); MLX contributors; HIPS Autograd contributors. Upstream model-card and code-license notices are under `work/sleep_frailty_ml/hyenadna/`. See the report for data and code licensing distinctions. No AlphaGenome outputs were used for training.
'''
paths=[]
for path in WORK.iterdir():
    if path.is_file() and path.suffix in ['.py','.md','.R']:paths.append(path)
for name in ['project_sequence_pairs.csv']:
    paths.append(WORK/name)
for name in ['model.safetensors','config.json','README.md','UPSTREAM_CODE_LICENSE.txt']:
    paths.append(WORK/'hyenadna'/name)
paths.extend(p for p in OUT.rglob('*') if p.is_file() and p.suffix!='.zip' and p.name not in ['package_receipt.json','finalization.log'])
paths.append(ROOT/'outputs/ncds_score_variant_request.csv')
digest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
archive=OUT/'sleep_frailty_ai_v1.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    z.writestr('sleep_frailty_ai_v1/START_HERE.md',start)
    z.writestr('sleep_frailty_ai_v1/package_checksums.json',json.dumps(digest,indent=2)+'\n')
    for path in sorted(set(paths)):z.write(path,'sleep_frailty_ai_v1/'+str(path.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
(OUT/'package_receipt.json').write_text(json.dumps({'archive':archive.name,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'file_count':len(set(paths))+2},indent=2)+'\n')
print('Verified model package',archive,archive.stat().st_size,'bytes')
