import csv,hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent;previous=p.parents[1]/'physical_function_ldsc/insomnia_source_sensitivity/prepare_receipt.json';prior=json.loads(previous.read_text());n=prior['common_valid_nonpalindromic'];fall_receipt=json.loads((p.parent/'download_receipt.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''):h.update(c)
 return h.hexdigest()
assert sha(fall_receipt['path'])==fall_receipt['sha256']
for key in ('questionnaire','registry'):assert sha(prior['output_paths'][key])==prior['output_sha256'][key]
values={};allrg={};reported={}
for name in ('questionnaire','registry'):
 t=(p/f'{name}.log').read_text();assert 'Analysis finished at' in t and int(re.search(r'(\d+) SNPs with valid alleles',t).group(1))==n
 rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',t).groups());pv=float(re.search(r'Z-score: [\d.-]+\s+P: ([\d.eE-]+)',t).group(1))
 prefix=p/f'{name}falls.ldsc.bgz_{name}.sumstats.gz';v={k:[float(i) for i in Path(str(prefix)+'.'+k+'.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')};assert all(len(a)==200 for a in v.values()) and all(x>0 and math.isfinite(x) for k in ('hsq1','hsq2') for x in v[k])
 d=[z/math.sqrt(a*b) for z,a,b in zip(v['gencov'],v['hsq1'],v['hsq2'])];m=sum(d)/200;sd=math.sqrt(199/200*sum((x-m)**2 for x in d));assert abs(sd-se)<.00005
 values[name]=v;allrg[name]=d;reported[name]={'rg':rg,'se':se,'p':pv}
assert values['questionnaire']['hsq1']==values['registry']['hsq1']
diff=reported['questionnaire']['rg']-reported['registry']['rg'];ds=[x-y for x,y in zip(allrg['questionnaire'],allrg['registry'])];av=sum(ds)/200;se=math.sqrt(199/200*sum((x-av)**2 for x in ds));pv=math.erfc(abs(diff/se)/math.sqrt(2));out={'common_markers':n,'individual':reported,'contrast':{'questionnaire_minus_registry':diff,'paired_se':se,'p':pv,'ci95':[diff-1.95996398454*se,diff+1.95996398454*se]},'blocks':200,'fall_heritability_delete_vectors_identical':True,'questionnaire_registry_input_sha256_verified':True,'note':'Contrast uses rounded log rg; a post-selected sensitivity, not a causal comparison.'}
(p/'results.json').write_text(json.dumps(out,indent=2)+'\n')
with (p/'rg_deletes.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['block','questionnaire','registry']);w.writerows((i+1,allrg['questionnaire'][i],allrg['registry'][i]) for i in range(200))
print(json.dumps(out,indent=2))
