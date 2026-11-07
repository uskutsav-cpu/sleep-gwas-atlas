import hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent
out={}
for trait in ('questionnaire','registry','walking','grip'):
 receipt=json.loads((p/f'stream_receipt_{trait}.json').read_text());assert receipt['source_md5']=='d9b5c1e7513e209bc57326264c95c356' and receipt['source_bytes']==269614908 and receipt['gzip_members']==282 and receipt['retained']==984470 and receipt['duplicate_eligible_rows']==0
 t=(p/f'run2_{trait}.log').read_text();assert 'Analysis finished at' in t and int(re.search(r'(\d+) SNPs with valid alleles',t).group(1))==984470
 rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',t).groups());pvalue=float(re.search(r'Z-score: [\d.-]+\s+P: ([\d.eE-]+)',t).group(1))
 block=t.split('Genetic Covariance')[1].split('Genetic Correlation')[0];i,ise=map(float,re.search(r'Intercept: ([\d.-]+) \(([\d.]+)\)',block).groups())
 out[trait]={'rg':rg,'se':se,'p_logged':pvalue,'p_four_test_bonferroni':None if pvalue==0 else min(1,4*pvalue),'p_zero_is_numerical_underflow':pvalue==0,'cross_trait_intercept':i,'cross_trait_intercept_se':ise,'log_sha256':hashlib.sha256(t.encode()).hexdigest(),'n_snps':984470}
(p/'summary.json').write_text(json.dumps({'traits':out,'warning':'Feasibility only: factor inputs, factor direction and N semantics are unverified; no validated sleep-independent frailty GWAS.'},indent=2)+'\n');print(json.dumps(out,indent=2))
