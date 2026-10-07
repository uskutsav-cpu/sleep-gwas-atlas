#!/usr/bin/env python3
"""Rebuild exact float64 signedLD from small versioned dosage factors; no network."""
import argparse,csv,hashlib,json,pathlib
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent
RAW=HERE.parents[3]/'work'/'ld_genotypes_research_v1'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for d in iter(lambda:f.read(1024*1024),b''):h.update(d)
 return h.hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=pathlib.Path,default=RAW);args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
 results=[]
 for c in [5,6,11]:
  diag=json.loads((HERE/f'chr{c}_reconstruction_diagnostics.json').read_text())['summary'];factor=HERE/f'chr{c}_genotype_factors.npz';assert sha(factor)==diag['genotype_factor_sha256'];f=np.load(factor);rows=list(csv.DictReader((HERE/f'chr{c}_variants.tsv').open(),delimiter='\t'));keys=np.array([r['variant_key'] for r in rows]);assert np.array_equal(f['variant_keys'],keys);G=f['ALT_dosage'];centered=G.astype(np.float64)-G.mean(axis=0);X=centered/np.sqrt(np.sum(centered**2,axis=0));R=X.T@X
  temp=args.output_dir/f'chr{c}_rebuild_validation.npz';np.savez_compressed(temp,R=R,variant_keys=keys,effect_allele=np.array([r['ALT'] for r in rows]));digest=sha(temp);assert digest==diag['matrix_sha256'],f'Exact matrix rebuild hash mismatch{c}'
  target=args.output_dir/f'chr{c}_signed_LD.npz'
  if target.exists():assert sha(target)==digest;temp.unlink()
  else:temp.rename(target)
  result={'chr':c,'samples':G.shape[0],'SNPs':G.shape[1],'expected_matrix_sha256':diag['matrix_sha256'],'recomputed_matrix_sha256':digest,'status':'PASS_EXACT_FLOAT64_FACTOR_REBUILD','no_original_raw_VCF_required':True};results.append(result);print(json.dumps(result),flush=True)
 (HERE/'factor_rebuild_validation.json').write_text(json.dumps({'script_sha256':sha(pathlib.Path(__file__)),'numpy':np.__version__,'results':results},indent=2)+'\n')
if __name__=='__main__':main()
