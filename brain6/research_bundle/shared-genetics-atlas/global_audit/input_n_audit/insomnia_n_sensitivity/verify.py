import hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent
files=['build.py','build_receipt.json','run.sh','rg_constantN.log','rg_rawN.log','constantN_console.log','rawN_console.log','matched_constantN.sumstats.gz','matched_rawN.sumstats.gz']
def sha(x):
 h=hashlib.sha256()
 with x.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
results={}
for name in ('constantN','rawN'):
 s=(p/f'rg_{name}.log').read_text();assert 'Analysis finished at' in s
 marker=re.search(r'(\d+) SNPs with valid alleles',s);rg=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',s);assert marker and rg
 results[name]={'markers':int(marker.group(1)),'rg':float(rg.group(1)),'se':float(rg.group(2)),'z':float(rg.group(3)),'p':float(rg.group(4))}
assert results['constantN']['markers']==results['rawN']['markers']==1054978
assert results['constantN']['rg']==.4033 and results['constantN']['se']==.0318 and results['constantN']['p']==5.836e-37
out={'results':results,'files_sha256':{x:sha(p/x) for x in files},'notes':'Numerical check only, same-data post-hoc N sensitivity, not proof either N convention is appropriate.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(results,indent=2))
