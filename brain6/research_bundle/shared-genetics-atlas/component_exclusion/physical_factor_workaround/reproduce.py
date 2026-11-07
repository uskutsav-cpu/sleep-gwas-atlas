"""Reproduce saved-summary analyses in a new directory; no SSD or downloads."""
import argparse,hashlib,json,shutil,subprocess,sys
from pathlib import Path

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--output-dir',required=True,type=Path,help='New directory whose parent already exists; existing paths are refused.')
 args=ap.parse_args();src=Path(__file__).resolve().parent;out=args.output_dir.resolve()
 if out.exists():raise SystemExit('Refusing existing output path: '+str(out))
 if not out.parent.is_dir():raise SystemExit('Output parent must already exist')
 rscript=shutil.which('Rscript')
 if not rscript:raise SystemExit('Rscript is required; base R is sufficient')
 import numpy  # Confirm the only Python dependency before creating output.
 pairs=[('walking','grip'),('walking','activity'),('grip','activity')]+[(t,s) for s in ('insomnia','questionnaire') for t in ('walking','grip','activity')]
 scripts=['joint_prepare.py','joint_fit.R','residual_fit.R','identification_range.py','explicit_index.py','weight_sensitivity.py','draw_index.R']
 names=scripts+['factor_feasibility_insomnia.json','factor_feasibility_questionnaire.json']
 for a,b in pairs:
  name=a+'_'+b;names.append(name+'.log')
  names.extend(name+a+'.common.sumstats.gz_'+b+'.common.sumstats.gz.'+key+'.delete' for key in ('hsq1','hsq2','gencov'))
 missing=[name for name in names if not (src/name).is_file()]
 if missing:raise SystemExit('Missing required files: '+', '.join(missing))
 source_hashes={name:sha(src/name) for name in names}
 out.mkdir()
 for name in names:shutil.copy2(src/name,out/name)
 steps=[('prepare',[sys.executable,'joint_prepare.py']),('joint_fit',[rscript,'joint_fit.R']),('residual_fit',[rscript,'residual_fit.R']),('identification',[sys.executable,'identification_range.py']),('index',[sys.executable,'explicit_index.py']),('weights',[sys.executable,'weight_sensitivity.py']),('figure',[rscript,'draw_index.R'])]
 for name,cmd in steps:
  with (out/(name+'.stdout.log')).open('w') as stdout,(out/(name+'.stderr.log')).open('w') as stderr:
   subprocess.run(cmd,cwd=out,stdout=stdout,stderr=stderr,check=True)
 matches={}
 for name in ('explicit_index.json','identification_ranges.json','weight_sensitivity.json'):
  matches[name]=json.loads((src/'joint_fit'/name).read_text())==json.loads((out/'joint_fit'/name).read_text())
 if not all(matches.values()):raise RuntimeError('Replayed numerical JSON differs from saved results: '+str(matches))
 if any(sha(src/name)!=value for name,value in source_hashes.items()):raise RuntimeError('A source input changed during replay')
 result={'source_directory':str(src),'output_directory':str(out),'source_sha256':source_hashes,'exact_numerical_json_matches':matches,'python':sys.version,'numpy':numpy.__version__,'rscript':rscript,'scope':'Saved-summary calculations and figures only; original GWAS processing, LDSC fits, phenotype metadata and independent validation are not reproduced.'}
 (out/'replay_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
