import hashlib,json,re
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; out=p/'joint_fit';out.mkdir(exist_ok=True)
pairs=[('walking','grip'),('walking','activity'),('grip','activity')]+[(t,s) for s in ('insomnia','questionnaire') for t in ('walking','grip','activity')]
known={}
for s in ('insomnia','questionnaire'):
 d=json.loads((p/('factor_feasibility_'+s+'.json')).read_text())
 for r in d['pairs']:known[r['pair']]=r['log_sha256']
values=[];deletes=[];traits={};hashes={}
for a,b in pairs:
 name=a+'_'+b;log=p/(name+'.log');raw=log.read_bytes();hashes[log.name]=hashlib.sha256(raw).hexdigest();assert hashes[log.name]==known[name]
 text=raw.decode();assert 'Analysis finished at' in text
 assert int(re.search(r'(\d+) SNPs with valid alleles',text).group(1))==1053614
 point,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
 base=p/(name+a+'.common.sumstats.gz_'+b+'.common.sumstats.gz')
 arrays={}
 for key in ('hsq1','hsq2','gencov'):
  path=Path(str(base)+'.'+key+'.delete');data=path.read_bytes();hashes[path.name]=hashlib.sha256(data).hexdigest();arrays[key]=np.array([float(x) for x in data.split()]);assert len(arrays[key])==200 and np.isfinite(arrays[key]).all()
 for trait,key in ((a,'hsq1'),(b,'hsq2')):
  assert (arrays[key]>0).all()
  if trait in traits:assert np.array_equal(traits[trait],arrays[key])
  else:traits[trait]=arrays[key]
 d=arrays['gencov']/np.sqrt(arrays['hsq1']*arrays['hsq2']);assert abs(np.sqrt(199/200*np.sum((d-d.mean())**2))-se)<.00005
 values.append(point);deletes.append(d)
deletes=np.asarray(deletes).T;centered=deletes-deletes.mean(axis=0);cov=199/200*(centered.T@centered);eig=np.linalg.eigvalsh(cov);assert eig.min()>0
np.savetxt(out/'moments.txt',values);np.savetxt(out/'deletes.txt',deletes);np.savetxt(out/'covariance.txt',cov)
receipt={'pairs':[a+'_'+b for a,b in pairs],'source_hashes':hashes,'minimum_sampling_eigenvalue':float(eig.min()),'sampling_condition_number':float(eig.max()/eig.min()),'n_blocks':200,'snp_count':1053614,'per_trait_heritability_vectors_match':True}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_hashes'},indent=2))
