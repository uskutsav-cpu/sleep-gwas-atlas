"""Finish checks, reports and packaging after a specific live training process."""
from pathlib import Path
import sys,os,time,argparse,subprocess,json,hashlib,csv
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent;OUT=ROOT/'outputs/sleep_frailty_ai_v1'
p=argparse.ArgumentParser();p.add_argument('--wait-for-pid',type=int);args=p.parse_args()
if args.wait_for_pid:
    print('Waiting for the verified training PID',args.wait_for_pid,flush=True)
    while True:
        ps=subprocess.run(['ps','-p',str(args.wait_for_pid),'-o','command='],capture_output=True,text=True)
        if ps.returncode!=0 or 'sleep_frailty_ml/train.py' not in ps.stdout:break
        time.sleep(10)
    print('Training process is terminal; checking saved results',flush=True)
if not (OUT/'test_results.json').is_file():raise RuntimeError('Training stopped without final test results; inspect training.log')
def run(script,python=sys.executable):
    print('RUN',script,flush=True);subprocess.run([python,str(WORK/script)],cwd=ROOT,check=True)
for script in ['merge_cnn_results.py','select_model.py','evaluate_uncertainty.py','check_inference.py']:run(script)
selected=json.loads((OUT/'model_selection.json').read_text())['selected_model']
checkpoints={'kmer_logistic':'kmer_baseline.npz','frozen_hyenadna':'frozen_hyenadna.npz','finetuned_hyenadna':'finetuned_hyenadna.npz','small_dna_cnn':'cnn.npz'}
for stem,checkpoint in [('project_candidate_predictions','finetuned_hyenadna.npz'),('project_cnn_predictions','cnn.npz'),('project_recommended_model_predictions',checkpoints[selected])]:
    target=OUT/(stem+'.csv');weights=OUT/'models'/checkpoint
    if target.exists():
        receipt=json.loads(target.with_suffix('.provenance.json').read_text())
        assert receipt['checkpoint_sha256']==hashlib.sha256(weights.read_bytes()).hexdigest()
        assert receipt['output_sha256']==hashlib.sha256(target.read_bytes()).hexdigest()
        print('Verified existing predictions',target.name,flush=True)
    else:subprocess.run([sys.executable,str(WORK/'predict.py'),'--pairs',str(WORK/'project_sequence_pairs.csv'),'--checkpoint',str(weights),'--output',str(target),'--backend','mlx'],cwd=ROOT,check=True)
run('export_embeddings.py')
subprocess.run(['Rscript',str(WORK/'independent_metrics.R'),str(OUT)],check=True)
with (OUT/'independent_R_metrics.csv').open() as f:independent=list(csv.DictReader(f))
results=json.loads((OUT/'test_results.json').read_text())['models']
differences={}
for row in independent:
    differences[row['model']]={}
    for metric in ['auroc','average_precision']:
        delta=abs(float(row[metric])-results[row['model']][metric]);differences[row['model']][metric]=delta
        # Training-side float32 integration and R float64 rank arithmetic differ at ~1e-9.
        assert delta<1e-7,(row['model'],metric,delta)
(OUT/'independent_metric_agreement.json').write_text(json.dumps({'absolute_differences':differences,'tolerance':1e-7,'reason':'Float32 training metric integration versus independent double-precision R rank arithmetic'},indent=2)+'\n')
run('build_report.py')
plot_python=Path('/Users/swethasunilkumar/Documents/Codex/2026-10-29/investigate-whether-we-can-directly-reconstruct-2/work/venv/bin/python3')
if plot_python.is_file():run('figures.py',str(plot_python))
else:print('Plotting environment absent; run figures.py with matplotlib installed',flush=True)
(OUT/'COMPLETED.json').write_text(json.dumps({'status':'completed experiment, not clinically validated','selected_model':selected,'independent_R_metrics_agree':True,'model_package':'sleep_frailty_ai_v1.zip'},indent=2)+'\n')
run('package.py')
print('FINALIZATION COMPLETE',flush=True)
