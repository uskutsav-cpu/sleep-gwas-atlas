import csv,collections,math,pathlib,json,hashlib
p=pathlib.Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(p/'predictions.tsv');b=read(p.parent/'masked_predictions.tsv');s=read(p/'summary.tsv');t=read(p/'training_objectives.tsv');w=read(p/'fold_weights.tsv');assert len(a)==len(b)==1080 and len(t)==250 and len(w)==50
key=lambda x:(x['repeat'],x['fold'],x['sleep_trait'],x['outcome']);prior={key(x):x for x in b};assert len(prior)==1080 and len({key(x) for x in a})==1080
for x in a:
 z=prior[key(x)];assert x['observed']==z['observed'] and float(x['unweighted_rank2'])==float(z['rank2_multistart'])
for group in set((x['repeat'],x['fold']) for x in t):
 vals=[x for x in t if (x['repeat'],x['fold'])==group]; chosen=[x for x in vals if x['selected']=='True'];assert len(chosen)==1 and float(chosen[0]['objective'])==min(float(x['objective']) for x in vals)
for x in s:
 rows=[v for v in a if x['group']=='ALL' or v['outcome']==x['group']];e=[float(v[x['model']])-float(v['observed']) for v in rows];se=[float(v['SE']) for v in rows];assert len(rows)==int(x['n']);assert math.isclose(math.sqrt(sum(z*z for z in e)/len(e)),float(x['rmse']),rel_tol=1e-12);assert math.isclose(math.sqrt(sum((z/q)**2 for z,q in zip(e,se))/len(e)),float(x['se_scaled_rmse']),rel_tol=1e-12)
result={'status':'PASS','same_masks_and_baseline_predictions':True,'best_training_objectives':True,'summary_metrics_reproduced':True,'files_sha256':{name:hashlib.sha256((p/name).read_bytes()).hexdigest() for name in ['fit_weighted.py','predictions.tsv','training_objectives.tsv','fold_weights.tsv','summary.tsv']},'caveat':'Arithmetic and mask checks do not prove calibration or independent validation.'};(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS same masks, training selection, pooled and per-outcome metrics')
