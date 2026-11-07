import csv,hashlib,json,math,re
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;atlas=p.parents[1];base=atlas/'multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run';old=base/'jackknife_uncertainty';prev=list(old.glob('ldsc-*.log'));new=list(p.glob('ldsc-*.log'));assert len(prev)==len(new)==1
ss=[x.read_text() for x in (prev[0],new[0])];assert all('Analysis finished' in x for x in ss) and ss[0].count('991707 SNPs with valid alleles.')==9 and ss[1].count('991707 SNPs with valid alleles.')==3
pat=r'Genetic Correlation: ([\d.\-]+) \(([\d.]+)\)\nZ-score: [\d.\-]+\nP: ([\deE.\-]+)';a,b=[re.findall(pat,x) for x in ss];assert len(a)==9 and len(b)==3
names=['FI','FFS','HFRS'];points_old=dict(zip(['adhd','mdd','scz','bipolar','alz','parkinson']+names,[(float(x[0]),float(x[1])) for x in a]));points_new=dict(zip(names,[(float(x[0]),float(x[1])) for x in b]));files=[prev[0],new[0],base/'input_verification.json'];r=json.loads((base/'input_verification.json').read_text());assert len(r['files'])==10 and all(z['rows']==991707 for z in r['files'].values())
def vec(directory,first,second,pt):
 parts=[]
 for typ in ('hsq1','hsq2','gencov'):
  matches=list(directory.glob('*'+first+'.sumstats.gz_'+second+'.sumstats.gz.'+typ+'.delete'));assert len(matches)==1;files.append(matches[0]);v=np.loadtxt(matches[0]);assert v.shape==(200,) and np.isfinite(v).all();parts.append(v)
 h1,h2,cov=parts;assert (h1>0).all() and (h2>0).all();rg=cov/np.sqrt(h1*h2);se=float(np.sqrt(199/200*np.sum((rg-rg.mean())**2)));assert abs(se-pt[1])<0.00006;return rg,h1,h2
ia,ih,ah=vec(old,'MVP','adhd',points_old['adhd']);out=[];min_eigen=1.
for n in names:
 sf,ih2,fh=vec(old,'MVP',n,points_old[n]);af,ah2,fh2=vec(p,'adhd',n,points_new[n]);assert np.array_equal(ih,ih2) and np.array_equal(ah,ah2) and np.array_equal(fh,fh2)
 S,A,F=points_old['adhd'][0],points_new[n][0],points_old[n][0]
 def calc(x,y,z):
  d=(1-x*x)*(1-y*y);assert np.all(d>0);v=(z-x*y)/np.sqrt(d);assert np.all(np.abs(v)<=1);return v
 point=float(calc(S,A,F));deletes=calc(ia,af,sf)
 for x,y,z in [(S,A,F),*zip(ia,af,sf)]:
  e=float(np.linalg.eigvalsh([[1,x,z],[x,1,y],[z,y,1]]).min());min_eigen=min(min_eigen,e);assert e>=-1e-10
 se=float(np.sqrt(199/200*np.sum((deletes-deletes.mean())**2)));pv=math.erfc(abs(point/se)/math.sqrt(2));out.append({'frailty':n,'insomnia_adhd_rg':S,'adhd_frailty_rg':A,'insomnia_frailty_rg':F,'partial_genetic_rg':point,'paired_jackknife_se':se,'normal_95_low':point-1.96*se,'normal_95_high':point+1.96*se,'unadjusted_p':pv,'bonferroni_p_three':min(1,3*pv)})
with (p/'results.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'matched_snps':991707,'blocks':200,'smallest_three_trait_matrix_eigenvalue':min_eigen,'identical_delete_vectors':'Insomnia, ADHD, and each frailty trait have identical corresponding h2 deletion vectors across the two runs','input_sha256':{str(x.relative_to(atlas)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files},'qualification':'conditional genetic correlation, not causal interaction or individual-level conditioning'},indent=2)+'\n');print(json.dumps(out,indent=2));print('minimum eigenvalue',min_eigen)
