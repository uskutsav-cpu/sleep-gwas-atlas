import csv,hashlib,json,math,re
from pathlib import Path
p=Path(__file__).resolve().parent
prep=json.loads((p/'prepare_receipt.json').read_text()); assert prep['ordered_keys_check']
assert json.loads((p/'duplicate_check.json').read_text())['duplicated_ids']==0
for f,expected in prep['output_sha256'].items():
 h=hashlib.sha256()
 with open(f,'rb') as x:
  for b in iter(lambda:x.read(1024*1024),b''):h.update(b)
 assert h.hexdigest()==expected
vals={};deleted={};pts={};ses={}
for t in ('FI','FFS'):
 text=(p/(t+'.log')).read_text(); assert 'Analysis finished at' in text
 assert int(re.search(r'After merging with summary statistics, (\d+) SNPs remain',text).group(1))==prep['stats']['matched']
 pts[t],loggedse=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',text).groups())
 base=p/(t+t+'.sumstats.gz_insomnia.sumstats.gz')
 vals[t]={k:[float(x) for x in Path(str(base)+'.'+k+'.delete').read_text().split()] for k in ('hsq1','hsq2','gencov')}
 assert all(len(v)==200 for v in vals[t].values())
 assert all(math.isfinite(x) and x>0 for k in ('hsq1','hsq2') for x in vals[t][k])
 d=[c/math.sqrt(a*b) for a,b,c in zip(vals[t]['hsq1'],vals[t]['hsq2'],vals[t]['gencov'])]
 assert all(math.isfinite(x) for x in d)
 avg=sum(d)/200;se=math.sqrt(199/200*sum((x-avg)**2 for x in d));assert abs(se-loggedse)<0.00005
 deleted[t]=d;ses[t]=se
assert vals['FI']['hsq2']==vals['FFS']['hsq2']
d=[a-b for a,b in zip(deleted['FI'],deleted['FFS'])];avg=sum(d)/200;se=math.sqrt(199/200*sum((v-avg)**2 for v in d));delta=pts['FI']-pts['FFS'];z=delta/se;pval=math.erfc(abs(z)/math.sqrt(2))
res={'n_snps':prep['stats']['matched'],'fi_rg_printed':pts['FI'],'ffs_rg_printed':pts['FFS'],'fi_jackknife_se':ses['FI'],'ffs_jackknife_se':ses['FFS'],'contrast':'FI minus FFS','difference_rounded_points':delta,'paired_se':se,'ci95_unadjusted_low':delta-1.95996398454*se,'ci95_unadjusted_high':delta+1.95996398454*se,'z':z,'p_unadjusted':pval,'blocks':200,'insomnia_delete_vectors_identical':True,'limitations':'Point difference uses rounded r_g values; post-hoc single contrast is not adjusted for the atlas-wide family. Shared support and input identity checks do not establish cohort identity or correct trait model.'}
(p/'paired_result.json').write_text(json.dumps(res,indent=2)+'\n')
with (p/'delete_rg.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['block','FI_rg_delete','FFS_rg_delete','difference']);w.writerows((i+1,deleted['FI'][i],deleted['FFS'][i],d[i]) for i in range(200))
print(json.dumps(res,indent=2))
