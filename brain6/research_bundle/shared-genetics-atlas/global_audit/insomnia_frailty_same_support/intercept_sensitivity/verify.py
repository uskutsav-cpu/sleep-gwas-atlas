import hashlib,json,math,pathlib,re
p=pathlib.Path(__file__).resolve().parent;rows=[]
for fixed in ('0.165','0.193','0'):
 f=p/f'intercept_{fixed}.log';text=f.read_text();assert '1134953 SNPs with valid alleles.' in text and 'Analysis finished at' in text and f'Intercept: constrained to {fixed}' in text
 m=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',text)
 if m:
  rg,se,z,pvalue=map(float,m.groups());rows.append({'fixed_intercept':fixed,'status':'ESTIMATED','rg':rg,'se':se,'p':pvalue})
 else:
  assert 'rg out of bounds' in text and fixed=='0';rows.append({'fixed_intercept':fixed,'status':'RG_OUT_OF_BOUNDS','rg':None,'se':None,'p':None})
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
r={'rows':rows,'common_snps':1134953,'logs_sha256':{x:sha(p/f'intercept_{x}.log') for x in ('0.165','0.193','0')},'script_sha256':sha(p/'run.sh'),'limitation':'Fixed-intercept stress is not an estimator of the true intercept, calibrated uncertainty, independent replication or evidence against cohort overlap.'};(p/'receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(rows,indent=2))
