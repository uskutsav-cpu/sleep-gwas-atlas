#!/usr/bin/env python3
import csv,gzip,math,json,pathlib,hashlib,collections
from contextlib import ExitStack
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation')
files={n:base/'common_support'/f'{n}.sumstats.gz' for n in ('S','FI','FFS')};files['HFRS']=base/'filtered/HFRS_HM3.tsv';files['position_FI']=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz')
out=pathlib.Path('outputs/shared-genetics-atlas/four_trait_marginal');n=1138262;cut=.05/(4*n)
hfrs={};dups=set()
with files['HFRS'].open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  k=r['SNP']
  if k in hfrs:dups.add(k)
  else:hfrs[k]=(r['CHR'],r['A1'],r['A2'],r['BETA'],r['P'])
for k in dups:hfrs.pop(k,None)
counts=collections.Counter();seen=set();hits=[]
with ExitStack() as stack:
 fs=[stack.enter_context(gzip.open(files[x],'rt')) for x in ('S','FI','FFS')]
 readers=[csv.DictReader(f,delimiter='\t') for f in fs]
 for a,b,c in zip(*readers,strict=True):
  key=tuple(a[x] for x in ('SNP','A1','A2'));assert all(tuple(z[x] for x in ('SNP','A1','A2'))==key for z in (b,c));assert key[0] not in seen;seen.add(key[0]);z=[float(x['Z']) for x in (a,b,c)];assert all(math.isfinite(x) for x in z)
  ps=[math.erfc(max(abs(x)-.0005,0)/math.sqrt(2)) for x in z];passes=[name for name,p in zip(('insomnia','FI','FFS'),ps) if p<=cut]
  h=hfrs.get(key[0]);status='MISSING';aligned='';hp='';hb='';ha='';hother='';hchr=''
  if key[0] in dups:status='DUPLICATE'
  elif h:
   hchr,ha,hother,hb,hp=h
   if len({key[1],key[2]})!=2 or any(len(x)!=1 or x not in 'ACGT' for x in (key[1],key[2],ha,hother)) or {key[1],key[2]}!={ha,hother}:status='ALLELE_MISMATCH'
   else:
    try:
     hp=float(hp);hb=float(hb);assert math.isfinite(hp) and 0<hp<=1 and math.isfinite(hb)
     aligned=hb if (key[1],key[2])==(ha,hother) else -hb;status='MATCHED'
     if hp<=cut:passes.append('HFRS')
    except (ValueError,AssertionError):status='INVALID'
  counts['HFRS_'+status]+=1
  if passes:
   hits.append(dict(SNP=key[0],A1=key[1],A2=key[2],Z_insomnia=z[0],Z_FI=z[1],Z_FFS=z[2],p_insomnia_upper=ps[0],p_FI_upper=ps[1],p_FFS_upper=ps[2],HFRS_status=status,HFRS_CHR=hchr,HFRS_A1=ha,HFRS_A2=hother,HFRS_beta_A1=hb,HFRS_beta_common_A1=aligned,HFRS_p=hp,threshold_crossing='+'.join(passes),chr37='',bp37='',position_status=''))
assert len(seen)==n
wanted={r['SNP'] for r in hits};pos=collections.defaultdict(list)
with gzip.open(files['position_FI'],'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in wanted:pos[r['SNP']].append((r['CHR'],r['BP'],r['A1'],r['A2']))
for x in hits:
 match=pos[x['SNP']]
 if len(match)==1 and {match[0][2],match[0][3]}=={x['A1'],x['A2']}:
  x['chr37'],x['bp37']=match[0][:2];x['position_status']='UNIQUE_ALLELE_MATCH'
  if x['HFRS_status']=='MATCHED' and x['chr37']!=x['HFRS_CHR']:
   x['HFRS_status']='CHR_MISMATCH_AFTER_POSITION_CHECK';x['HFRS_beta_common_A1']='';x['threshold_crossing']='+'.join(t for t in x['threshold_crossing'].split('+') if t!='HFRS') or 'NONE_AFTER_HOLD'
 else:x['position_status']='MISSING_OR_AMBIGUOUS'
 counts['pattern_'+x['threshold_crossing']]+=1
 counts['selected_position_'+x['position_status']]+=1
 counts['selected_final_HFRS_'+x['HFRS_status']]+=1
with (out/'crossings.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(hits[0]),delimiter='\t');w.writeheader();w.writerows(hits)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while v:=f.read(4*1024*1024):h.update(v)
 return h.hexdigest()
summary={'n_support':n,'family_size':4*n,'threshold':cut,'counts':dict(counts),'selected_count':len(hits),'inputs':{k:{'path':str(v),'sha256':sha(v)} for k,v in files.items()},'caveat':'p-based marginal patterns only; no independent loci, causal effects, specificity or independent replication; HFRS not queried outside 3-trait support.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
