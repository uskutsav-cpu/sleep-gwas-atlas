#!/usr/bin/env python3
"""Implementation arithmetic checks only, not phenotype validation."""
from pathlib import Path
import csv,hashlib,json,subprocess,sys,tempfile
import numpy as np
from score_15_research import calculate,read_rows,WEIGHTS
O=Path(__file__).resolve().parent;W=O.parent/'work';model=read_rows(O/'unvalidated_15_variant_score_prototype.csv')
with np.load(W/'stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True) as f:
    G=np.asarray(f['genotypes'],dtype=float);ids=f['variants'].astype(str)
assert list(ids)==[r['variant'] for r in model] and G.shape==(503,15) and np.isfinite(G).all()
obs=[]
for i,g in enumerate(G):
 for j,r in enumerate(model):
  flip=j%2==1;obs.append({'sample_id':f'{i:04d}','variant':r['variant'],'dosage_allele':r['other_allele'] if flip else r['counted_allele'],'other_allele':r['counted_allele'] if flip else r['other_allele'],'dosage':2-g[j] if flip else g[j]})
res=calculate(model,obs,allow_verified_reversal=True);actual=np.array([[r['FI_raw_sum'],r['FFS_raw_sum_inferred'],r['insomnia_raw_sum']] for r in res]);expected=G@np.array([[float(r[k]) for k in WEIGHTS] for r in model]);maxerr=float(np.max(np.abs(expected-actual)));assert maxerr<1e-12
checks={}
for name,rows in [('missing',obs[:-1]),('duplicate',obs+[obs[0]]),('empty',[]),('bad_dosage',[{**obs[0],'dosage':float('nan')}]+obs[1:]),('mismatched_allele',[{**obs[0],'dosage_allele':'N'}]+obs[1:])]:
 try:calculate(model,rows,allow_verified_reversal=True)
 except ValueError:checks[name]=True
 else:raise AssertionError(f'{name} incorrectly accepted')
with tempfile.TemporaryDirectory() as tmp:
 p=Path(tmp);inp=p/'bad.csv';out=p/'out.csv'
 with inp.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['sample_id','variant','dosage_allele','other_allele','dosage']);w.writeheader();w.writerow(obs[0])
 cmd=[sys.executable,str(O/'score_15_research.py'),'--dosages',str(inp),'--output',str(out),'--allow-verified-reversal'];r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode!=0 and not out.exists();checks['cli_missing_no_output']=True
 with inp.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['sample_id','variant','dosage_allele','other_allele','dosage']);w.writeheader();w.writerows(obs[:15])
 r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode==0 and out.exists();out.write_text('preserve');r=subprocess.run(cmd,capture_output=True,text=True);assert r.returncode!=0 and out.read_text()=='preserve';checks['valid_cli_write_and_existing_output_preserved']=True
r={'reference':'1000 Genomes EUR 503, genotypes already aligned to counted alleles','variants':15,'samples':503,'input_contrasts_deliberately_reversed':7,'max_absolute_arithmetic_difference':maxerr,'negative_checks':checks,'scorer_sha256':hashlib.sha256((O/'score_15_research.py').read_bytes()).hexdigest(),'weights_sha256':hashlib.sha256((O/'unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'limitation':'Implementation consistency only: no participant outcome, reference or study allele validity, target sample independence, score performance, or clinical validation.'};(O/'score_15_reference_implementation_check.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
