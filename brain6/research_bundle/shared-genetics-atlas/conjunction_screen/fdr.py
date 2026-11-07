import csv,gzip,math,json,pathlib,bisect,hashlib
root=pathlib.Path('outputs/shared-genetics-atlas/conjunction_screen'); base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support'); names=['S','FI','FFS'];fs=[base/(x+'.sumstats.gz') for x in names]; ps=[];seen=set()
with gzip.open(fs[0],'rt') as a,gzip.open(fs[1],'rt') as b,gzip.open(fs[2],'rt') as c:
 for rows in zip(*(csv.DictReader(f,delimiter='\t') for f in (a,b,c)),strict=True):
  s=rows[0];key=tuple(s[k] for k in ('SNP','A1','A2'));assert all(tuple(t[k] for k in ('SNP','A1','A2'))==key for t in rows);assert key[0] not in seen;seen.add(key[0])
  z=[float(t['Z']) for t in rows];assert all(math.isfinite(x) for x in z); pv=[math.erfc(max(abs(x)-.0005,0)/math.sqrt(2)) for x in z];ps.extend((max(pv[0],pv[1]),max(pv[0],pv[2])))
m=len(ps);assert m==2*1138262;ps.sort();harmonic=sum(1/i for i in range(1,m+1));bh_count=by_count=0
for j,p in enumerate(ps,1):
 if p<=.05*j/m:bh_count=j
 if p<=.05*j/(m*harmonic):by_count=j
# BH adjusted p at selected p values: q_i = min_{j>=i} m*p_j/j, clipped at 1.
top=list(csv.DictReader((root/'top50_per_pair.tsv').open(),delimiter='\t'));indices={bisect.bisect_left(ps,float(r['p_conjunction'])) for r in top};q={};qmin=1.0
for i in range(m-1,-1,-1):
 qmin=min(qmin,m*ps[i]/(i+1))
 if i in indices:q[i]=qmin
for r in top:
 i=bisect.bisect_left(ps,float(r['p_conjunction']));assert ps[i]==float(r['p_conjunction']);r['q_BH']=q[i];r['q_BY']=min(1,q[i]*harmonic)
with (root/'top50_fdr.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(top[0]),delimiter='\t');w.writeheader();w.writerows(top)
summary={'number_hypotheses':m,'level':.05,'harmonic_number':harmonic,'BH_rejections':bh_count,'BY_rejections':by_count,'minimum_p':ps[0],'minimum_q_BH':qmin,'minimum_q_BY':min(1,qmin*harmonic),'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in fs}}
(root/'fdr_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
