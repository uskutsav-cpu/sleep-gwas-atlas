import json,math,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;folder=p/'joint_fit'
y=np.loadtxt(folder/'moments.txt');d=np.loadtxt(folder/'deletes.txt');assert y.shape==(9,) and d.shape==(200,9)
def index(row):
 variance=3+2*sum(row[:3]);assert variance>0
 for offset in (3,6):
  m=np.eye(4)
  for (a,b),v in zip(((0,1),(0,2),(1,2),(0,3),(1,3),(2,3)),np.r_[row[:3],row[offset:offset+3]]):m[a,b]=m[b,a]=v
  assert np.linalg.eigvalsh(m).min()>-1e-10
 return np.array([sum(row[3:6]),sum(row[6:9])])/math.sqrt(variance)
point=index(y);deletes=np.array([index(row) for row in d]);assert (abs(deletes)<=1).all() and (abs(point)<=1).all()
cov=199/200*((deletes-deletes.mean(axis=0)).T@(deletes-deletes.mean(axis=0)))
rows=[]
for i,source in enumerate(('registry','questionnaire')):
 se=math.sqrt(cov[i,i]);pv=math.erfc(abs(point[i]/se)/math.sqrt(2))
 rows.append({'source':source,'genetic_index_correlation':float(point[i]),'block_jackknife_se':se,'unadjusted_95_ci':[float(point[i]-1.95996398454*se),float(point[i]+1.95996398454*se)],'normal_p':pv,'bonferroni_two_sources':min(1,2*pv)})
delta=point[1]-point[0];dv=deletes[:,1]-deletes[:,0];se=math.sqrt(199/200*np.sum((dv-dv.mean())**2))
result={'weights':[1,1,1],'index_genetic_variance_before_normalization':float(3+2*sum(y[:3])),'estimates':rows,'paired_questionnaire_minus_registry':{'estimate':float(delta),'block_jackknife_se':se,'normal_p':math.erfc(abs(delta/se)/math.sqrt(2))},'n_blocks':200,'moments_sha256':hashlib.sha256((folder/'moments.txt').read_bytes()).hexdigest(),'delete_values_sha256':hashlib.sha256((folder/'deletes.txt').read_bytes()).hexdigest(),'qualification':'Post-hoc equal-weight index of standardized genetic components, not a measured physical-frailty score or matched exhaustion-deletion GWAS. No zero indicator-residual covariance restriction. Source/model assumptions remain.'}
(folder/'explicit_index.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
