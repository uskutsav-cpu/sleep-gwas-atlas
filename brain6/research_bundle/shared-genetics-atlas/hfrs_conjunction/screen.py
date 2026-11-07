#!/usr/bin/env python3
import csv,hashlib,json,math,heapq,collections,array
from pathlib import Path
base=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered')
out=Path('outputs/shared-genetics-atlas/hfrs_conjunction')
paths={'insomnia':base/'R13_INSOMNIA_HM3.tsv','hfrs':base/'HFRS_HM3.tsv'}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
def parse(p):
 with p.open() as f:
  yield from csv.DictReader(f,delimiter='\t')
# One slot per source row remains in the family, even if excluded.
s={};duplicates=set();n=0
for r in parse(paths['insomnia']):
 n+=1;k=r['SNP']
 if k in s:duplicates.add(k)
 else:s[k]=r
for k in duplicates:s.pop(k,None)
# Count all HFRS matches first, ensuring ambiguous source duplicates are excluded.
counts=collections.Counter()
for r in parse(paths['hfrs']):
 if r['SNP'] in s:counts[r['SNP']]+=1
status=collections.Counter();status['duplicate_insomnia_unique_ids']=len(duplicates)
ps=array.array('d');top=[];hits=[];serial=0
for b in parse(paths['hfrs']):
 k=b['SNP'];a=s.get(k)
 if a is None:continue
 if counts[k]!=1:continue
 if a['CHR']!=b['CHR']:
  status['chromosome_mismatch']+=1;continue
 aa=(a['A1'].upper(),a['A2'].upper());bb=(b['A1'].upper(),b['A2'].upper())
 if len(set(aa))!=2 or any(x not in 'ACGT' or len(x)!=1 for x in aa+bb) or set(aa)!=set(bb):
  status['allele_mismatch_or_non_snp']+=1;continue
 try:
  pa,pb=float(a['P']),float(b['P']); ba,bbeta=float(a['BETA']),float(b['BETA'])
  if not all(math.isfinite(v) for v in (pa,pb,ba,bbeta)) or not (0<pa<=1 and 0<pb<=1):raise ValueError()
 except (ValueError,TypeError):status['invalid_p_or_beta']+=1;continue
 p=max(pa,pb);ps.append(p);serial+=1
 aligned=bbeta if aa==bb else -bbeta
 item={'snp':k,'chromosome':a['CHR'],'insomnia_A1':aa[0],'insomnia_A2':aa[1],'hfrs_A1':bb[0],'hfrs_A2':bb[1],'insomnia_beta_A1':ba,'hfrs_beta_A1':bbeta,'hfrs_beta_aligned_insomnia_A1':aligned,'insomnia_p':pa,'hfrs_p':pb,'conjunction_max_p':p,'bonferroni_p':min(1,p*n)}
 entry=(-p,serial,item)
 if len(top)<100:heapq.heappush(top,entry)
 elif p < -top[0][0]:heapq.heapreplace(top,entry)
 if p<=0.05/n:hits.append(item)
status['eligible']=len(ps);status['missing_hfrs_unique_insomnia']=sum(1 for k in s if counts[k]==0);status['duplicate_hfrs_unique_insomnia']=sum(1 for k in s if counts[k]>1)
ps=sorted(ps);harmonic=sum(1/i for i in range(1,n+1))
def nreject(c):
 j=0
 for i,p in enumerate(ps,1):
  if p<=.05*i/(n*c):j=i
 return j
fields=['snp','chromosome','insomnia_A1','insomnia_A2','hfrs_A1','hfrs_A2','insomnia_beta_A1','hfrs_beta_A1','hfrs_beta_aligned_insomnia_A1','insomnia_p','hfrs_p','conjunction_max_p','bonferroni_p']
for filename,items in [('top100.tsv',[x[2] for x in sorted(top,reverse=True)]),('hits.tsv',sorted(hits,key=lambda x:x['conjunction_max_p']))]:
 with (out/filename).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(items)
summary={'source_files':{k:{'path':str(v),'sha256':digest(v),'bytes':v.stat().st_size} for k,v in paths.items()},'family_size_insomnia_rows':n,'bonferroni_threshold':.05/n,'status':dict(status),'bonferroni_crossing_rows':len(hits),'BH_rejections_padded_family':nreject(1),'BY_rejections_padded_family':nreject(harmonic),'min_conjunction_p':ps[0] if ps else None,'note':'Exploratory SNP-level conjunction; source p-values as displayed; not independent cohorts, PLACO, conjFDR or colocalization.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
