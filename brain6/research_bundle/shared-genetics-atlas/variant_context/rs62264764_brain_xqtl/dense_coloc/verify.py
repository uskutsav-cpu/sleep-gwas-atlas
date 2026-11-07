import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;receipt=json.loads((p/'matching_receipt.json').read_text())
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as x:
  for b in iter(lambda:x.read(2**20),b''):h.update(b)
 return h.hexdigest()
def lse(x):
 m=max(x);return m+math.log(sum(math.exp(z-m) for z in x))
def labf(b,se,prior):
 v=se*se;r=prior**2/(prior**2+v);return .5*(math.log1p(-r)+r*(b/se)**2)
results=list(csv.DictReader((p/'posteriors.tsv').open(),delimiter='\t'));mx=0
for r in results:
 pair=r['cell']+'_'+r['trait']; info=receipt['pairs'][pair]; assert sha(info['file'])==info['sha256']
 d=list(csv.DictReader(open(info['file']),delimiter='\t'));assert len(d)==int(r['n'])>=500
 x=[labf(float(t['qtl_beta']),float(t['qtl_se']),.15) for t in d];y=[labf(float(t['gwas_beta_alt']),float(t['gwas_se']),.2) for t in d];xy=[a+b for a,b in zip(x,y)]
 a=lse(x);b=lse(y);c=lse(xy);ab=a+b
 h3=ab+math.log1p(-math.exp(c-ab))+2*math.log(1e-4)
 logs=[0,math.log(1e-4)+a,math.log(1e-4)+b,h3,math.log(float(r['p12']))+c];den=lse(logs);v=[math.exp(z-den) for z in logs]
 for k in range(5):mx=max(mx,abs(v[k]-float(r['H'+str(k)])))
assert len(results)==4 and mx<1e-10
old=receipt['qtl_source']['file_md5'];pub=list(csv.DictReader((p.parents[1]/'chr5_brain_xqtl/dense_coloc/publisher_snuc_manifest.tsv').open(),delimiter='\t'))
for kind in ('oli',):
 for suffix in ('.bed.gz','.bed.gz.tbi'):
  match=[r for r in pub if 'ROSMAP_CUIMC1_2_MIT' in r['File name'] and ('.'+kind+'.') in r['File name'] and r['File name'].endswith('_3_snp_all.20260330'+suffix)]
  assert len(match)==1 and match[0]['md5Sum']==old[kind+suffix]
out={'pairs':len(results),'max_absolute_posterior_difference':mx,'publisher_manifest_checks':2,'inputs_sha256_verified':True,'interpretation':'Arithmetic and data-integrity checks do not validate model assumptions or establish causality.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
