import csv,hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent;old=p.parent;prev=json.loads((old/'paired_receipt.json').read_text());scores={};v={};results=[]
for t in ('FI','FFS','HFRS'):
 text=(p/(t+'.log')).read_text();assert 'Analysis finished at' in text;assert int(re.search(r'After merging with summary statistics, (\d+) SNPs remain',text).group(1))==1099360
 point,se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups());pv=float(re.search(r'^P: ([\deE+.-]+)',text.split('Genetic Correlation\n')[1],re.M).group(1));base=p/(t+t+'.sumstats.gz_insomnia.sumstats.gz')
 vals={k:[float(x) for x in Path(str(base)+'.'+k+'.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')};assert all(len(x)==200 for x in vals.values());assert all(x>0 for k in ('hsq1','hsq2') for x in vals[k]);v[t]=vals
 d=[c/math.sqrt(a*b) for a,b,c in zip(vals['hsq1'],vals['hsq2'],vals['gencov'])];avg=sum(d)/200;se_j=math.sqrt(199/200*sum((x-avg)**2 for x in d));assert abs(se_j-se)<.00005;scores[t]=(point,d)
 oldrow=next(r for r in prev['correlations'] if r['trait']==t);results.append({'trait':t,'old_rg':oldrow['rg'],'new_rg':point,'old_se':oldrow['se_jackknife'],'new_se':se_j,'new_p':pv,'new_bonferroni_p_3':min(1,3*pv)})
assert v['FI']['hsq2']==v['FFS']['hsq2']==v['HFRS']['hsq2']
contrasts=[]
for a,b in [('HFRS','FI'),('HFRS','FFS'),('FI','FFS')]:
 d=[x-y for x,y in zip(scores[a][1],scores[b][1])];avg=sum(d)/200;se=math.sqrt(199/200*sum((x-avg)**2 for x in d));delta=scores[a][0]-scores[b][0];pv=math.erfc(abs(delta/se)/math.sqrt(2));oldrow=next(x for x in prev['contrasts'] if x['contrast']==a+' minus '+b)
 contrasts.append({'contrast':a+' minus '+b,'old_difference':oldrow['difference'],'new_difference':delta,'new_se_paired':se,'new_p_bonferroni_3':min(1,3*pv)})
for name,data in [('correlations.tsv',results),('contrasts.tsv',contrasts)]:
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
r={'results':results,'contrasts':contrasts,'markers':1099360,'num_delete_blocks':200,'identical_insomnia_delete_values':True,'prior_receipt_sha256':hashlib.sha256((old/'paired_receipt.json').read_bytes()).hexdigest(),'warning':'Post-hoc same-data N sensitivity; supplied Ns not validated.'}
(p/'paired_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
