import csv,gzip,math,hashlib,json,pathlib
base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support');dest=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/outputs/shared-genetics-atlas/placo_diagnostic'); dest.mkdir(exist_ok=True,parents=True)
files={k:base/(k+'.sumstats.gz') for k in ('S','FI','FFS')};counts={'n':0,'extreme':{'FI':0,'FFS':0},'zero':{'FI':0,'FFS':0}}; accum={k:{'v':[],'c':[]} for k in ('FI','FFS')}; sample=[];seen=set()
with gzip.open(files['S'],'rt') as a,gzip.open(files['FI'],'rt') as b,gzip.open(files['FFS'],'rt') as c:
 rs=[csv.DictReader(x,delimiter='\t') for x in (a,b,c)]
 for rows in zip(*rs,strict=True):
  s,fi,ffs=rows; key=(s['SNP'],s['A1'],s['A2']);assert all((x['SNP'],x['A1'],x['A2'])==key for x in rows); assert s['SNP'] not in seen;seen.add(s['SNP']); counts['n']+=1
  z=[float(x['Z']) for x in rows];assert all(math.isfinite(x) for x in z); ps=[math.erfc(abs(x)/math.sqrt(2)) for x in z]
  for idx,k in [(1,'FI'),(2,'FFS')]:
   if z[0]**2>80 or z[idx]**2>80:counts['extreme'][k]+=1;continue
   if z[0]==0 or z[idx]==0:counts['zero'][k]+=1
   if not (ps[0]<1e-4 and ps[idx]<1e-4):accum[k]['v'].append((z[0],z[idx]))
   if ps[0]>=1e-4 and ps[idx]>=1e-4:accum[k]['c'].append((z[0],z[idx]))
  if int.from_bytes(hashlib.blake2b(s['SNP'].encode(),digest_size=8,person=b'placoqc1').digest(),'big')%512==0:
   sample.append([s['SNP'],s['A1'],s['A2'],*z])
with (dest/'random_controls.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['SNP','A1','A2','Z_insomnia','Z_FI','Z_FFS']);w.writerows(sample)
def moments(a):
 n=len(a);mx=sum(x for x,y in a)/n;my=sum(y for x,y in a)/n;xx=sum((x-mx)**2 for x,y in a)/(n-1);yy=sum((y-my)**2 for x,y in a)/(n-1);xy=sum((x-mx)*(y-my) for x,y in a)/(n-1);return {'n':n,'var1':xx,'var2':yy,'cor':xy/math.sqrt(xx*yy)}
params={k:{'variance_subset':moments(d['v']),'correlation_subset':moments(d['c'])} for k,d in accum.items()}
record={'source':'R13 FinnGen insomnia vs original FI/FFS, HM3 common LDSC support, NOT dense whole-genome variants','input_hashes':{k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in files.items()},'counts':counts,'parameters':params,'sample_n':len(sample),'sample_rule':'BLAKE2b SNP id hash mod 512=0; independent of effect and p','p_from_z':'two sided normal tail, rounded Z to 0.001 from LDSC munged inputs','z2_exclusion':'either z squared > 80','parameter_threshold':1e-4}
(dest/'diagnostic_parameters.json').write_text(json.dumps(record,indent=2)+'\n'); print(json.dumps(record,indent=2))
