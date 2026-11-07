import csv,hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent; earlier=p.parent/'paired_common_support'; match=json.loads((p/'matching_receipt.json').read_text()); assert match['ordered_four_way_keys_identical']
file=Path('/Volumes/Brain6Work/brain6-frailty-paired/HFRS.sumstats.gz')
h=hashlib.sha256()
with file.open('rb') as f:
 for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
assert h.hexdigest()==match['hfrs_input_sha256']
vals={};deletes={};point={};ses={};ps={}
for t in ('FI','FFS','HFRS'):
 folder=p if t=='HFRS' else earlier
 text=(folder/(t+'.log')).read_text();assert 'Analysis finished at' in text
 assert int(re.search(r'After merging with summary statistics, (\d+) SNPs remain',text).group(1))==match['matching_rows']
 point[t],se=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
 ps[t]=float(re.search(r'^P: ([\deE+.-]+)',text.split('Genetic Correlation\n')[1],re.M).group(1))
 base=folder/(t+t+'.sumstats.gz_insomnia.sumstats.gz')
 vals[t]={k:[float(x) for x in Path(str(base)+'.'+k+'.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')}
 assert all(len(z)==200 for z in vals[t].values());assert all(math.isfinite(x) and x>0 for k in ('hsq1','hsq2') for x in vals[t][k])
 d=[c/math.sqrt(a*b) for a,b,c in zip(vals[t]['hsq1'],vals[t]['hsq2'],vals[t]['gencov'])];assert all(math.isfinite(x) for x in d)
 m=sum(d)/200;ses[t]=math.sqrt(199/200*sum((v-m)**2 for v in d));assert abs(ses[t]-se)<.00005;deletes[t]=d
assert vals['FI']['hsq2']==vals['FFS']['hsq2']==vals['HFRS']['hsq2']
contrasts=[]
for a,b in [('HFRS','FI'),('HFRS','FFS'),('FI','FFS')]:
 d=[x-y for x,y in zip(deletes[a],deletes[b])];m=sum(d)/200;se=math.sqrt(199/200*sum((v-m)**2 for v in d));delta=point[a]-point[b];z=delta/se;pv=math.erfc(abs(z)/math.sqrt(2))
 contrasts.append({'contrast':a+' minus '+b,'difference':delta,'se_paired':se,'ci95_unadjusted_low':delta-1.95996398454*se,'ci95_unadjusted_high':delta+1.95996398454*se,'z':z,'p_unadjusted':pv,'bonferroni_p_3':min(1,3*pv)})
results=[{'trait':t,'rg':point[t],'se_jackknife':ses[t],'p':ps[t],'bonferroni_p_3':min(1,3*ps[t])} for t in ('FI','FFS','HFRS')]
for name,rows in [('contrasts.tsv',contrasts),('correlations.tsv',results)]:
 with (p/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
with (p/'delete_rg.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['block','FI','FFS','HFRS']);w.writerows([i+1,deletes['FI'][i],deletes['FFS'][i],deletes['HFRS'][i]] for i in range(200))
receipt={'n_snps':match['matching_rows'],'all_insomnia_hsq_delete_vectors_identical':True,'num_blocks':200,'correlations':results,'contrasts':contrasts,'p_family':'3 post-hoc correlations and 3 post-hoc contrasts each corrected separately; not atlas-wide','limitation':'Point differences use four-decimal log r_g values; results subject to support, phenotype, ancestry, overlap and source-model assumptions.'}
(p/'paired_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
