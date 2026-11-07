import csv,json,hashlib,pathlib,numpy as np
p=pathlib.Path(__file__).resolve().parent;source=p.parent/'ukb_cross_tab.csv';fp=['robust','prefrail','frail'];fi=['robust','mild','moderate_severe'];matrix=np.zeros((3,3),dtype=np.int64);seen=set()
with source.open() as f:
 for r in csv.DictReader(f):
  i=fp.index(r['fried_category']);j=fi.index(r['fi_category']);assert (i,j) not in seen;seen.add((i,j));matrix[i,j]=int(r['n'])
assert len(seen)==9 and matrix.sum()==383845 and np.all(matrix>=0)
dist=np.abs(np.arange(3)[:,None]-np.arange(3)[None,:]);weights={'unweighted':(dist==0).astype(float),'linear':1-dist/2,'quadratic':1-(dist/2)**2}
def measures(m):
 q=m/m.sum(); marg=q.sum(axis=1)[:,None]*q.sum(axis=0)[None,:];o={}
 for name,w in weights.items():
  obs=float((w*q).sum());exp=float((w*marg).sum());o[name]=(obs-exp)/(1-exp);o[name+'_observed']=obs;o[name+'_expected']=exp
 return o
observed=measures(matrix);rng=np.random.default_rng(20261027); draws=rng.multinomial(int(matrix.sum()),(matrix/matrix.sum()).ravel(),size=10000);boot={name:np.array([measures(a.reshape(3,3))[name] for a in draws]) for name in weights}
results=[]
for name in weights:
 lo,hi=np.quantile(boot[name],[.025,.975]);results.append({'weight':name,'n':int(matrix.sum()),'observed_agreement':observed[name+'_observed'],'expected_agreement':observed[name+'_expected'],'kappa':observed[name],'bootstrap_2.5pct':float(lo),'bootstrap_97.5pct':float(hi),'bootstrap_n':len(draws)})
with (p/'agreement.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
receipt={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256((p/'calculate.py').read_bytes()).hexdigest(),'n':int(matrix.sum()),'bootstrap_seed':20261027,'bootstrap_replicates':len(draws),'rng':'numpy default_rng PCG64','numpy':np.__version__,'weight_definitions':{'unweighted':'1 if equal else 0','linear':'1 - |i-j|/2','quadratic':'1 - ((i-j)/2)^2'},'caveat':'Categories differ; weights arbitrary; independent multinomial bootstrap is descriptive'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(results)
