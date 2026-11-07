import csv,gzip,math,json,pathlib,hashlib,heapq
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support')
out=pathlib.Path('outputs/shared-genetics-atlas/conjunction_screen');n=1138262; threshold=.05/(2*n)
paths=[base/(x+'.sumstats.gz') for x in ('S','FI','FFS')]; counts={'FI':0,'FFS':0}; hits=[];best={'FI':[],'FFS':[]};seen=set();total=0
with gzip.open(paths[0],'rt') as a,gzip.open(paths[1],'rt') as b,gzip.open(paths[2],'rt') as c:
 for s,fi,ffs in zip(*(csv.DictReader(f,delimiter='\t') for f in (a,b,c)),strict=True):
  key=tuple(s[k] for k in ('SNP','A1','A2'));assert all(tuple(r[k] for k in ('SNP','A1','A2'))==key for r in (fi,ffs));assert s['SNP'] not in seen;seen.add(s['SNP']);total+=1
  z=[float(r['Z']) for r in (s,fi,ffs)];assert all(math.isfinite(x) for x in z)
  ps=[math.erfc(max(abs(x)-.0005,0)/math.sqrt(2)) for x in z]
  for i,pair in enumerate(('FI','FFS'),1):
   p=max(ps[0],ps[i])
   row=[pair,*key,z[0],z[i],ps[0],ps[i],p,min(1,p*2*n)]
   entry=(-p,total,row)
   if len(best[pair])<50:heapq.heappush(best[pair],entry)
   elif p < -best[pair][0][0]:heapq.heapreplace(best[pair],entry)
   if p<=threshold:counts[pair]+=1;hits.append(row)
assert total==n
with (out/'hits.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['pair','SNP','A1','A2','Z_insomnia','Z_frailty','p_insomnia_upper','p_frailty_upper','p_conjunction','p_bonferroni_2pairs']);w.writerows(sorted(hits,key=lambda r:r[-2]))
with (out/'top50_per_pair.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['pair','SNP','A1','A2','Z_insomnia','Z_frailty','p_insomnia_upper','p_frailty_upper','p_conjunction','p_bonferroni_2pairs']);w.writerows(sorted([item[2] for heap in best.values() for item in heap],key=lambda r:(r[0],r[-2])))
meta={'n_markers':n,'number_of_hypotheses':2*n,'threshold':threshold,'rounding_margin':0.0005,'counts':counts,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'limitations':'post-selected exploratory common HapMap3 only; rounded p-derived Z; no independent locus count'}
(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
