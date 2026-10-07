#!/usr/bin/env python3
"""Independent Python coloc posterior calculation from native SuSiE Bayes factors."""
import csv,hashlib,json,math,pathlib
import numpy as np
from scipy.special import logsumexp
HERE=pathlib.Path(__file__).resolve().parent;LOCAL=HERE.parents[3]/'work/ld_genotypes_research_v1/C_finemap_inputs'
def weights(a,b,p1,p2,p12):
 l1=float(logsumexp(a));l2=float(logsumexp(b));l4=float(logsumexp(a+b));l3=l1+l2+math.log1p(-math.exp(l4-l1-l2))
 logs=np.array([0,math.log(p1)+l1,math.log(p2)+l2,math.log(p1)+math.log(p2)+l3,math.log(p12)+l4]);return np.exp(logs-logsumexp(logs))
def main():
 rows=list(csv.DictReader((HERE/'native_coloc_all_signal_pairs.tsv').open(),delimiter='\t'));results=[];sources={}
 for x in rows:
  path=LOCAL/(x['model']+'_native_BFs_LOCAL_ONLY.tsv');z=list(csv.DictReader(path.open(),delimiter='\t'));a=np.array([float(r['logBF1']) for r in z]);b=np.array([float(r['logBF2']) for r in z]);post=weights(a,b,float(x['p1']),float(x['p2']),float(x['p12']));native=np.array([float(x[f'PP.H{i}.abf']) for i in range(5)]);err=float(np.max(np.abs(post-native)));assert err<1e-12,(x['model'],err)
  results.append({'model':x['model'],'p12':x['p12'],'n_snps':len(z),'native_H4':native[4],'independent_python_H4':post[4],'max_posterior_abs_error_H0_H4':err,'status':'PASS_INDEPENDENT_NUMERICAL_RECOMPUTATION','scope':'METHOD_NUMERICAL_VALIDATION_NOT_INDEPENDENT_COHORT_OR_MODEL_CALIBRATION'})
  sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
 with (HERE/'independent_coloc_recomputation.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
 (HERE/'independent_coloc_validation_provenance.json').write_text(json.dumps({'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'native_table_sha256':hashlib.sha256((HERE/'native_coloc_all_signal_pairs.tsv').read_bytes()).hexdigest(),'input_bf_hashes':sources,'native_R_BFs_used':True,'independent_posterior_implementation':True,'number_recomputed':len(results),'max_error':max(r['max_posterior_abs_error_H0_H4'] for r in results)},indent=2)+'\n')
 print(json.dumps({'recomputed':len(results),'max_abs_error':max(r['max_posterior_abs_error_H0_H4'] for r in results)}))
if __name__=='__main__':main()
