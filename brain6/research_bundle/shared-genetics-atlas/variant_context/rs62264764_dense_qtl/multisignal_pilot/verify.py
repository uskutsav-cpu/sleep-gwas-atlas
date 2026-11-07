import pathlib,csv,json,math,hashlib,numpy as np,gzip
p=pathlib.Path(__file__).resolve().parent
with (p/'variants.tsv').open() as f:v=list(csv.DictReader(f,delimiter='\t'))
with gzip.open(p/'ld_matrix.tsv.gz','rt') as f:R=np.loadtxt(f,delimiter='\t')
assert R.shape==(len(v),len(v)) and len(v)==1548 and np.isfinite(R).all();assert np.max(np.abs(R-R.T))<1e-9 and np.max(np.abs(np.diag(R)-1))<1e-9
with (p/'single_signal_same_support.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
def lsum(x):
 m=max(x);return m+math.log(math.fsum(math.exp(y-m) for y in x))
def bf(b,se,sd):
 r=sd*sd/(sd*sd+se*se);return .5*(math.log1p(-r)+r*(b/se)**2)
maxerr=0
for x in rows:
 t=x['trait'];a=[bf(float(r['qtl_beta']),float(r['qtl_se']),.15) for r in v];b=[bf(float(r[t+'_beta']),float(r[t+'_se']),.2) for r in v];la=lsum(a);lb=lsum(b);lab=lsum([a0+b0 for a0,b0 in zip(a,b)]);ld=la+lb+math.log1p(-math.exp(lab-la-lb));h=[0,math.log(1e-4)+la,math.log(1e-4)+lb,math.log(1e-8)+ld,math.log(float(x['p12']))+lab];den=lsum(h);z=[math.exp(q-den) for q in h];err=max(abs(z[i]-float(x['H'+str(i)])) for i in range(5));assert err<1e-10;maxerr=max(maxerr,err)
assert len(rows)==4
with (p/'credible_sets.tsv').open() as f:cs=list(csv.DictReader(f,delimiter='\t'))
assert len(cs)==2 and all(len(x['variants'].split(','))==int(x['size']) for x in cs)
result={'n_variants':len(v),'matrix_max_asymmetry':float(np.max(np.abs(R-R.T))),'matrix_max_diagonal_error':float(np.max(np.abs(np.diag(R)-1))),'single_signal_rows_independently_verified':len(rows),'maximum_single_signal_error':maxerr,'credible_set_membership_counts_verified':True,'single_signal_sha256':hashlib.sha256((p/'single_signal_same_support.tsv').read_bytes()).hexdigest(),'note':'arithmetic and matrix-format checks, not LD consistency or statistical calibration'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
