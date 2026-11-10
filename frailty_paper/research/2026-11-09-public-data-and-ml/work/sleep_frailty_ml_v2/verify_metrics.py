from common import *
import subprocess

sets=[('reused_benchmark',OUT/'reused_benchmark_predictions.npz',OUT/'reused_benchmark_results.json'),('external',WORK/'external_predictions_local.npz',OUT/'external_assay_results.json')]
agreement={}
for name,path,result in sets:
    if not path.exists():continue
    saved=np.load(path);model_names=[k for k in saved.files if k not in ['label','row_index']];y=saved['label'];tmp=WORK/(name+'_metric_check_local.csv')
    with tmp.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['row_id','label',*model_names]);w.writerows([i,int(yy),*[float(saved[k][i]) for k in model_names]] for i,yy in enumerate(y))
    target=OUT/(name+'_independent_R_metrics.csv');subprocess.run(['Rscript',str(WORK/'independent_metrics.R'),str(tmp),str(target)],check=True)
    expected=json.loads(result.read_text())['models'];delta={}
    for row in csv.DictReader(target.open()):
        delta[row['model']]={key:abs(float(row[key])-expected[row['model']][key]) for key in ['auroc','average_precision']}
        assert max(delta[row['model']].values())<1e-10,(name,row['model'],delta[row['model']])
    agreement[name]=delta
(OUT/'independent_metric_agreement.json').write_text(json.dumps(agreement,indent=2)+'\n');print('Independent R metrics verified',list(agreement),flush=True)
