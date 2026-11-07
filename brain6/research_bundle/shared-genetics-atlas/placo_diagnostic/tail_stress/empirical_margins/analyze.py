#!/usr/bin/env python3
import collections,csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;parent=p.parents[1]
params=json.loads((parent/'diagnostic_parameters.json').read_text())
base=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support')
files={k:base/(k+'.sumstats.gz') for k in ('S','FI','FFS')}
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
for k,path in files.items():assert sha(path)==params['input_hashes'][k]
acc={(pair,trait):{'n':0,'sum':0.,'sum2':0.,'counts':collections.Counter()} for pair in ('FI','FFS') for trait in ('S',pair)}
seen=set();total=0
with gzip.open(files['S'],'rt') as a,gzip.open(files['FI'],'rt') as b,gzip.open(files['FFS'],'rt') as c:
 readers=[csv.DictReader(f,delimiter='\t') for f in (a,b,c)]
 for s,fi,ffs in zip(*readers,strict=True):
  key=(s['SNP'],s['A1'],s['A2']);assert all((x['SNP'],x['A1'],x['A2'])==key for x in (fi,ffs));assert s['SNP'] not in seen;seen.add(s['SNP']);total+=1
  z={'S':float(s['Z']),'FI':float(fi['Z']),'FFS':float(ffs['Z'])};assert all(math.isfinite(x) for x in z.values())
  for pair in ('FI','FFS'):
   u,v=z['S'],z[pair]
   if u*u>80 or v*v>80 or (math.erfc(abs(u)/math.sqrt(2))<1e-4 and math.erfc(abs(v)/math.sqrt(2))<1e-4):continue
   for trait in ('S',pair):
    d=acc[(pair,trait)];x=z[trait];d['n']+=1;d['sum']+=x;d['sum2']+=x*x
    for t in (3,4,5,6):d['counts'][t]+=abs(x)>t
assert total==params['counts']['n']
out=[];summary={}
for (pair,trait),d in acc.items():
 expected=params['parameters'][pair]['variance_subset'];assert d['n']==expected['n']
 mean=d['sum']/d['n'];var=(d['sum2']-d['n']*mean*mean)/(d['n']-1);fitted=expected['var1' if trait=='S' else 'var2'];assert math.isclose(var,fitted,rel_tol=1e-10)
 summary[pair+'_'+trait]={'n':d['n'],'mean':mean,'variance':var,'fitted_variance':fitted}
 for t in (3,4,5,6):
  count=d['counts'][t];prob=math.erfc(t/math.sqrt(2*fitted));exp=d['n']*prob
  out.append({'pair':pair,'trait':trait,'n':d['n'],'abs_z_threshold':t,'observed_count':count,'normal_reference_expected_count':exp,'observed_div_expected':count/exp if exp else None})
with (p/'marginal_tails.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'input_sha256':params['input_hashes'],'total_rows':total,'subset_summary':summary,'output_sha256':sha(p/'marginal_tails.tsv'),'limitations':'Empirical mixture includes nonnull variants, ascertainment exclusion and LD. Normal reference is descriptive, not a type-I error test or evidence that any candidate is false.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
