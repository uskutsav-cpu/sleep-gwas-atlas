import csv,json,pathlib,math,hashlib,re
p=pathlib.Path(__file__).resolve().parent
with (p/'results.tsv').open() as f:points={r['trait']:float(r['restricted_rg']) for r in csv.DictReader(f,delimiter='\t')}
with (p.parent/'paired_contrasts.tsv').open() as f:old={r['contrast']:r for r in csv.DictReader(f,delimiter='\t')}
all_arrays={};deletes={};ses={}
for t in ('FI','FFS','HFRS'):
 base=p/f'{t}_deletions{t}.sumstats.gz_S.sumstats.gz';log=(p/f'{t}_deletions.log').read_text();assert '1134416 SNPs with valid alleles.' in log and 'Analysis finished at' in log
 match=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)',log);assert match and float(match.group(1))==points[t]
 v={k:[float(x) for x in pathlib.Path(str(base)+f'.{k}.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')}
 assert all(len(z)==200 for z in v.values()) and all(x>0 and math.isfinite(x) for k in ('hsq1','hsq2') for x in v[k]);r=[c/math.sqrt(a*b) for a,b,c in zip(v['hsq1'],v['hsq2'],v['gencov'])];assert all(math.isfinite(x) for x in r);m=sum(r)/200;se=math.sqrt(199/200*sum((x-m)**2 for x in r));assert abs(se-float(match.group(2)))<.00005
 all_arrays[t]=v;deletes[t]=r;ses[t]=se
assert all_arrays['FI']['hsq2']==all_arrays['FFS']['hsq2']==all_arrays['HFRS']['hsq2']
rows=[]
for a,b in [('HFRS','FI'),('HFRS','FFS'),('FI','FFS')]:
 delta=points[a]-points[b];vec=[x-y for x,y in zip(deletes[a],deletes[b])];mean=sum(vec)/200;se=math.sqrt(199/200*sum((x-mean)**2 for x in vec));pval=math.erfc(abs(delta/se)/math.sqrt(2));name=f'{a} minus {b}';rows.append({'contrast':name,'full_difference':old[name]['difference'],'restricted_difference':delta,'full_paired_se':old[name]['se_paired'],'restricted_paired_se':se,'restricted_ci95_low':delta-1.95996398454*se,'restricted_ci95_high':delta+1.95996398454*se,'restricted_p':pval,'restricted_bonferroni_p_3':min(1,3*pval)})
with (p/'paired_contrasts.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
inputs=[p/f'{t}_deletions{t}.sumstats.gz_S.sumstats.gz.{k}.delete' for t in ('FI','FFS','HFRS') for k in ('hsq1','hsq2','gencov')]+[p/f'{t}_deletions.log' for t in ('FI','FFS','HFRS')]+[p/'results.tsv',p.parent/'paired_contrasts.tsv',p/'support_receipt.json',p/'run_deletions.sh']
receipt={'blocks':200,'shared_insomnia_deletions_equal':True,'per_pair_se':ses,'input_sha256':{str(f):sha(f) for f in inputs},'caution':'Restricted point estimates use four-decimal logged rg; post-hoc sensitivity, not independent validation.'};(p/'paired_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(rows,indent=2))
