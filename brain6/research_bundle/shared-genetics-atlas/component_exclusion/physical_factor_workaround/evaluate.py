import hashlib,json,math,re,sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2

p=Path(__file__).resolve().parent
receipt=json.loads((p/'prepare_receipt.json').read_text()); data=Path(receipt['output_directory'])
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while True:
   b=f.read(1048576)
   if not b: break
   h.update(b)
 return h.hexdigest()
assert all(sha(data/(k+'.common.sumstats.gz'))==v for k,v in receipt['output_sha256'].items())
sleep=sys.argv[1] if len(sys.argv)>1 else 'insomnia'
assert sleep in ('insomnia','questionnaire')
if sleep=='questionnaire':
 qreceipt=json.loads((p/'questionnaire_receipt.json').read_text())
 assert sha(Path(qreceipt['output_path']))==qreceipt['output_sha256'] and qreceipt['selected_rows']==receipt['common_nonpalindromic_rows']
pairs=[('walking','grip'),('walking','activity'),('grip','activity'),('walking',sleep),('grip',sleep),('activity',sleep)]
points=[]; dels=[]; pertrait={}; summaries=[]
for a,b in pairs:
 name=a+'_'+b; t=(p/(name+'.log')).read_text()
 assert 'Analysis finished at' in t
 assert int(re.search(r'(\d+) SNPs with valid alleles',t).group(1))==receipt['common_nonpalindromic_rows']
 rg,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',t).groups())
 pv=float(re.search(r'Z-score: [\d.-]+\s+P: ([\d.eE-]+)',t).group(1))
 base=p/(name+a+'.common.sumstats.gz_'+b+'.common.sumstats.gz')
 v={k:np.loadtxt(str(base)+'.'+k+'.delete') for k in ['hsq1','hsq2','gencov']}
 assert all(len(x)==200 and np.isfinite(x).all() for x in v.values())
 assert np.all(v['hsq1']>0) and np.all(v['hsq2']>0)
 for trait,key in ((a,'hsq1'),(b,'hsq2')):
  if trait in pertrait: assert np.array_equal(pertrait[trait],v[key]),trait
  else: pertrait[trait]=v[key]
 d=v['gencov']/np.sqrt(v['hsq1']*v['hsq2']); jse=np.sqrt(199/200*np.sum((d-d.mean())**2))
 assert abs(jse-se)<.00005
 points.append(rg);dels.append(d);summaries.append({'pair':name,'rg':rg,'se':se,'p':pv,'log_sha256':sha(p/(name+'.log'))})
points=np.asarray(points); dels=np.asarray(dels).T
def loadings(x):
 wg,wa,ga=x[:3]
 if min(wg,wa,ga)<=0: return None
 return np.sqrt([wg*wa/ga,wg*ga/wa,wa*ga/wg])
lam=loadings(points); assert lam is not None
block_lam=[loadings(x) for x in dels]
invalid=sum(x is None or np.any(x>=1) for x in block_lam)
matrix=np.eye(4)
for (a,b),v in zip([(0,1),(0,2),(1,2),(0,3),(1,3),(2,3)],points): matrix[a,b]=matrix[b,a]=v
result={'pairs':summaries,'n_snps':receipt['common_nonpalindromic_rows'],'n_blocks':200,
        'heritability_delete_vectors_identical_per_trait':True,'loadings':lam.tolist(),
        'unique_variances':(1-lam**2).tolist(),'invalid_delete_loadings':int(invalid),
        'correlation_matrix_min_eigenvalue':float(np.linalg.eigvalsh(matrix).min()),
        'note':'Exploratory analytic one-factor check; paired LDSC blocks are not independent cohorts.'}
if np.any(lam>=1) or invalid or np.linalg.eigvalsh(matrix).min() < 0:
 result['gate']='failed: invalid factor loadings, leave-block fits, or observed correlation matrix'
else:
 q=points[3:]/lam
 dq=np.array([x[3:]/l for x,l in zip(dels,block_lam)])
 cov=199/200*((dq-dq.mean(axis=0)).T@(dq-dq.mean(axis=0)))
 one=np.ones(3); eig=np.linalg.eigvalsh(cov)
 if eig.min()<=0:
  result['gate']='failed: non-positive-definite projection covariance'
 else:
  inv=np.linalg.inv(cov); weights=inv@one/(one@inv@one); est=float(weights@q); diff=q-est
  Q=float(diff@inv@diff)
  result.update({'gate':'factor admissible under stated assumptions; heterogeneity test exploratory',
    'indicator_factor_correlations':q.tolist(),'projection_covariance':cov.tolist(),'gls_weights':weights.tolist(),
    'factor_insomnia_rg':est,'conditional_gls_se':float(np.sqrt(1/(one@inv@one))),
    'heterogeneity_chi2':Q,'heterogeneity_df':2,'heterogeneity_p':float(chi2.sf(Q,2)),
    'q_delete_min':dq.min(axis=0).tolist(),'q_delete_max':dq.max(axis=0).tolist()})
(p/('factor_feasibility_'+sleep+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2))
