import csv,json,pathlib,math,hashlib
p=pathlib.Path(__file__).resolve().parent
with (p/'correlations.tsv').open() as f:cor={x['frailty_trait']:x for x in csv.DictReader(f,delimiter='\t')}
old={t:json.loads((p.parents[1]/'frailty_sensitivity'/f'{t}_R13_ukbcommon.exact.json').read_text())[0] for t in ('FI','FFS')}
arrays={};ret={};se={}
for t in ('FI','FFS','HFRS'):
 base=p/f'{t}_deletions{t}.sumstats.gz_S.sumstats.gz'
 v={k:[float(x) for x in pathlib.Path(str(base)+f'.{k}.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')};assert all(len(z)==200 for z in v.values()) and all(x>0 and math.isfinite(x) for k in ('hsq1','hsq2') for x in v[k]);r=[c/math.sqrt(a*b) for a,b,c in zip(v['hsq1'],v['hsq2'],v['gencov'])];assert all(math.isfinite(x) for x in r)
 avg=sum(r)/200;se[t]=math.sqrt(199/200*sum((x-avg)**2 for x in r));assert abs(se[t]-float(cor[t]['se']))<.00005
 if t in old:
  o=[x[0] if isinstance(x,list) else x for x in old[t]['delete_rg']];assert len(o)==200 and max(abs(x-y) for x,y in zip(o,r))<1e-10
 arrays[t]=v;ret[t]=r
assert arrays['FI']['hsq2']==arrays['FFS']['hsq2']==arrays['HFRS']['hsq2']
points={t:(old[t]['rg'] if t in old else float(cor[t]['rg'])) for t in cor}
result=[]
for a,b in [('HFRS','FI'),('HFRS','FFS'),('FI','FFS')]:
 dif=[x-y for x,y in zip(ret[a],ret[b])];mean=sum(dif)/200;s=math.sqrt(199/200*sum((x-mean)**2 for x in dif));delta=points[a]-points[b];z=delta/s;pval=math.erfc(abs(z)/math.sqrt(2))
 result.append({'contrast':f'{a} minus {b}','difference':delta,'se_paired':s,'ci95_low':delta-1.95996398454*s,'ci95_high':delta+1.95996398454*s,'z':z,'p_two_sided':pval,'bonferroni_p_3':min(1,3*pval),'precision':'HFRS point estimate rounded to four decimals' if 'HFRS' in (a,b) else 'FI and FFS use stored exact points'})
with (p/'paired_contrasts.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
inputs=[p/f'{t}_deletions{t}.sumstats.gz_S.sumstats.gz.{k}.delete' for t in ('FI','FFS','HFRS') for k in ('hsq1','hsq2','gencov')]+[p/'correlations.tsv',p/'receipt.json',p/'run_deletions.sh']
receipt={'n_blocks':200,'shared_insomnia_delete_values_exactly_equal':True,'FI_and_FFS_match_prior_native_delete_rg':'within 1e-10','se':se,'p_family':'three post-hoc contrasts Bonferroni','source_sha256':{str(f):sha(f) for f in inputs},'limitation':'HFRS point rounded to 4 decimals; ancestry, cohort overlap, sample N and phenotype mismatch remain unresolved.'}
(p/'paired_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2))
