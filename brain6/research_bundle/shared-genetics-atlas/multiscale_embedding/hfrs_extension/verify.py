import csv,collections,math,pathlib,json,hashlib
p=pathlib.Path(__file__).resolve().parent
def read(n):
 with (p/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
r=read('masked_predictions.tsv');s=read('model_summary.tsv');t=read('training_objectives.tsv');src=read('source_matrix.tsv')
assert len(r)==1080 and len(src)==108 and len(t)==250
assert all(n==1 for n in collections.Counter((x['repeat'],x['sleep_trait'],x['outcome']) for x in r).values())
assert all(sum(x['selected']=='True' for x in t if (x['repeat'],x['fold'])==(rep,fold))==1 for rep in map(str,range(10)) for fold in map(str,range(5)))
for x in t:
 if x['selected']=='True':assert float(x['objective'])==min(float(z['objective']) for z in t if (z['repeat'],z['fold'])==(x['repeat'],x['fold']))
for row in s:
 z=[x for x in r if row['group']=='ALL' or x['outcome']==row['group']];v=math.sqrt(sum((float(x[row['model']])-float(x['observed']))**2 for x in z)/len(z));assert len(z)==int(row['n_predictions']) and math.isclose(v,float(row['rmse']),rel_tol=1e-12)
def sha(f):return hashlib.sha256((p/f).read_bytes()).hexdigest()
result={'status':'PASS','unique_masked_cell_per_repeat':True,'training_only_min_objective_checked':True,'all_summary_rmse_recomputed':True,'outputs_sha256':{x:sha(x) for x in ['source_matrix.tsv','masked_predictions.tsv','training_objectives.tsv','model_summary.tsv','coordinates.tsv','reconstruction.tsv','fit.py']},'caveat':'Checks result structure and arithmetic, not independent validation, input compatibility or causal interpretation.'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
