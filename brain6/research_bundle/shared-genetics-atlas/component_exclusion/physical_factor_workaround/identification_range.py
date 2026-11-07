import json,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;out={}
for source in ('insomnia','questionnaire'):
 d=json.loads((p/('factor_feasibility_'+source+'.json')).read_text());l=np.array(d['loadings']);u=1-l*l;r=np.array([x['rg'] for x in d['pairs'][3:]])
 assert (u>0).all()
 A=1+np.sum(l*l/u);B=np.sum(l*r/u);C=np.sum(r*r/u)-1
 discriminant=B*B-A*C;assert discriminant>=0
 low=(B-math.sqrt(discriminant))/A;high=(B+math.sqrt(discriminant))/A
 witnesses=[]
 candidates=[low,high]
 if low<=0<=high:candidates.append(0.)
 if high>0:candidates.append((max(0,low)+high)/2)
 if low<0:candidates.append((low+min(0,high))/2)
 for rho in candidates:
  c=(r-l*rho)/np.sqrt(u);remaining=1-rho*rho-c@c;assert remaining>-1e-10
  # Independent basis: F, three indicator-specific components, sleep-specific component.
  coefficients=np.zeros((5,5));coefficients[0,0]=1
  for i in range(3):coefficients[i+1,0]=l[i];coefficients[i+1,i+1]=np.sqrt(u[i])
  coefficients[4]=[rho,*c,math.sqrt(max(remaining,0))]
  gram=coefficients@coefficients.T
  assert np.max(abs(gram[1:4,4]-r))<1e-12
  assert np.max(abs(gram[0,1:4]-l))<1e-12
  assert np.min(np.linalg.eigvalsh(gram))>-1e-10
  assert np.max(abs(np.diag(gram)-1))<1e-10
  expected=np.eye(3)
  for a,b,v in ((0,1,d['pairs'][0]['rg']),(0,2,d['pairs'][1]['rg']),(1,2,d['pairs'][2]['rg'])):expected[a,b]=expected[b,a]=v
  assert np.max(abs(gram[1:4,1:4]-expected))<1e-12
  witnesses.append({'factor_sleep_correlation':float(rho),'indicator_residual_sleep_covariances':c.tolist(),'unallocated_sleep_variance':float(max(remaining,0)),'gram_min_eigenvalue':float(np.linalg.eigvalsh(gram).min()),'gram_matrix_F_walking_grip_activity_sleep':gram.tolist()})
 out[source]={'conditional_compatibility_range':[float(low),float(high)],'quadratic':[float(A),float(-2*B),float(C)],'loadings':l.tolist(),'sleep_indicator_correlations':r.tolist(),'witnesses':witnesses,'qualification':'Mathematical compatibility at rounded point moments under free indicator-specific residual covariance. Not a confidence interval, biological plausibility range, or estimate.'}
(p/'joint_fit/identification_ranges.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v['conditional_compatibility_range'] for k,v in out.items()},indent=2))
