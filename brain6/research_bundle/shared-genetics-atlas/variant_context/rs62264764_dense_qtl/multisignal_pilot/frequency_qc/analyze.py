import csv,pathlib,gzip,hashlib,math,json,collections,statistics
p=pathlib.Path(__file__).resolve().parent;pilot=p.parent
base=pathlib.Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data')
paths={t:base/'raw'/(t+'.txt.gz') for t in ('insomnia','adhd')}
expected={'adhd':'c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5','insomnia':'32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b'}
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as x:
  while b:=x.read(1048576):h.update(b)
 return h.hexdigest()
with (pilot/'variants.tsv').open() as f:variants=list(csv.DictReader(f,delimiter='\t'))
with (pilot/'sample_sizes.tsv').open() as f:ns={x['variant']:x for x in csv.DictReader(f,delimiter='\t')}
with (pilot/'kriging_diagnostics.tsv').open() as f:kr={(x['fit'],x['variant']):float(x['z_std_diff']) for x in csv.DictReader(f,delimiter='\t')}
assert len(variants)==1548 and len({x['rsid'] for x in variants})==1548
ids={x['rsid'] for x in variants};found={};source_rows={}
for trait in ('adhd','insomnia'):
 assert sha(paths[trait])==expected[trait]
 d=collections.defaultdict(list);n=0
 with gzip.open(paths[trait],'rt') as f:
  head=f.readline().split();ix={x:i for i,x in enumerate(head)}
  for line in f:
   n+=1;x=line.split();rs=x[ix['SNP']]
   if rs in ids:d[rs].append(dict(zip(head,x)))
 found[trait]=d;source_rows[trait]=n;print(trait,'source rows',n,'ids found',len(d),flush=True)
summary={};flags=[]
for trait in ('adhd','insomnia'):
 c=collections.Counter();diffs=[];groups=collections.defaultdict(list);n_deltas=[]
 for v in variants:
  rs=v['rsid'];matches=found[trait][rs]
  if len(matches)!=1:c['id_not_unique']+=1;continue
  x=matches[0];a,b=x['A1'].upper(),x['A2'].upper();ref,alt=v['ref'],v['alt']
  if x['CHR']!='3' or x['BP']!=v['pos37'] or {a,b}!={ref,alt}:c['coordinate_or_allele_mismatch']+=1;continue
  try:
   o=float(x['OR']);se=float(x['SE']);freq=float(x['FRQ_U_186843'] if trait=='adhd' else x['MAF']);beta=math.log(o)*(1 if a==alt else -1)
   if not (o>0 and se>0 and 0<=freq<=1 and all(map(math.isfinite,(beta,se,freq)))):raise ValueError()
  except (ValueError,KeyError):c['invalid_statistics']+=1;continue
  if not math.isclose(beta,float(v[trait+'_beta']),abs_tol=1e-10,rel_tol=1e-10) or not math.isclose(se,float(v[trait+'_se']),abs_tol=1e-10,rel_tol=1e-10):c['beta_or_se_mismatch']+=1;continue
  fref=float(v['reference_alt_frequency'])
  if trait=='adhd':source_alt=freq if a==alt else 1-freq;difference=abs(source_alt-fref)
  else:difference=abs(freq-min(fref,1-fref))
  z=kr[(trait,v['variant'])];diffs.append(difference);groups['kriging_abs_gt3' if abs(z)>3 else 'kriging_abs_le3'].append(difference);c['matched']+=1
  if difference>.05:c['difference_gt_005']+=1
  if difference>.10:c['difference_gt_010']+=1
  if difference>.15:c['difference_gt_015']+=1
  if trait=='adhd':
   ncase,ncontrol=float(x['Nca']),float(x['Nco']);neff=4*ncase*ncontrol/(ncase+ncontrol);n_deltas.append(abs(neff-float(ns[v['variant']]['adhd_N_AUDIT'])))
  if abs(z)>3 or difference>.10:flags.append({'trait':trait,'variant':v['variant'],'rsid':rs,'source_frequency_kind':'A1 control frequency aligned to ALT' if trait=='adhd' else 'source-labeled MAF','absolute_frequency_difference':difference,'kriging_z_std_diff':z})
 def stats(y):return {'n':len(y),'median_abs_diff':statistics.median(y) if y else None,'maximum_abs_diff':max(y) if y else None}
 summary[trait]={'counts':dict(c),'frequency_stats':stats(diffs),'groups':{k:stats(y) for k,y in groups.items()},'max_effective_n_difference':max(n_deltas) if n_deltas else None}
fields=['trait','variant','rsid','source_frequency_kind','absolute_frequency_difference','kriging_z_std_diff']
with (p/'flagged_variants.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(flags)
receipt={'sources':{t:{'path':str(paths[t]),'sha256':expected[t],'rows':source_rows[t]} for t in paths},'pilot_variants_sha256':sha(pilot/'variants.tsv'),'kriging_sha256':sha(pilot/'kriging_diagnostics.tsv'),'summary':summary,'flags_sha256':sha(p/'flagged_variants.tsv'),'caveat':'descriptive source and reference comparison, not evidence for or against allele coding, independent data or clinical use'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
