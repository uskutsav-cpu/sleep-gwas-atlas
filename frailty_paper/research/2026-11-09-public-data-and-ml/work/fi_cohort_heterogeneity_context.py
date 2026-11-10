import csv,gzip,hashlib,json,math,collections,pathlib,statistics
meta='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz';cohort='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz'
sample=lambda s:int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8,person=b'FIadapter-v1').digest(),'big')%64==0
m={};duplicates=set()
with gzip.open(meta,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if sample(r['SNP']):
   if r['SNP'] in m:duplicates.add(r['SNP'])
   else:m[r['SNP']]=r
for x in duplicates:m.pop(x,None)
records={};duplicates=set()
with gzip.open(cohort,'rt') as f:
 header=f.readline().rstrip('\r\n').split('\t'); idx=header.index('SNP')
 for line in f:
  fields=line.rstrip('\r\n').split('\t');rs=fields[idx];y=m.get(rs)
  if y is None:continue
  r=dict(zip(header,fields))
  if y is None or (r['CHR'],r['BP'])!=(y['CHR'],y['BP']) or {r['EA'],r['OA']}!={y['A1'],y['A2']} or r['EA']==r['OA']:continue
  if rs in records:duplicates.add(rs);continue
  try:
   u,su,t,st,b=map(float,(r['BETA_UKB1'],r['SE_UKB1'],r['BETA_TG'],r['SE_TG'],y['BETA']))
   if not all(map(math.isfinite,(u,su,t,st,b))) or min(su,st)<=0:continue
   sign=1 if r['EA']==y['A1'] else -1;u*=sign;t*=-sign
   records[rs]=(u,su,t,st,b,y['SE'])
  except ValueError:continue
for rs in duplicates:records.pop(rs,None)
prior=json.load(open('outputs/twingene_direction_genomewide_audit.json'))
assert len(records)==prior['counts']['valid_unique_matches'],(len(records),prior['counts']['valid_unique_matches'])
retained={};skips=collections.Counter()
for rs,(u,su,t,st,b,sm) in records.items():
 try:sm=float(sm)
 except ValueError:skips['bad_meta_se']+=1;continue
 if not math.isfinite(sm) or sm<=0:skips['bad_meta_se']+=1;continue
 weighted=(u/su**2+t/st**2)/(1/su**2+1/st**2)
 if abs(weighted-b)>0.0001:skips['reconstruction_mismatch']+=1;continue
 retained[rs]=((t-u)/math.hypot(su,st),b/sm,u,t)
records=retained
print('matched prior',prior['counts']['valid_unique_matches'],'eligible',len(records),'additional exclusions',dict(skips),flush=True)
ca=list(csv.DictReader(open('outputs/six_leads_fi_cohort_comparison.csv')));result={'sampling':'same deterministic 1-in-64 sample as twingene_direction_genomewide_audit','n_sample':len(records),'previous_adapter_audit_n':prior['counts']['valid_unique_matches'],'additional_exclusions':dict(skips),'groups':{},'candidates':{},'limitations':'Phenotypes and alleles differ and TwinGene sign remains inferred; SNPs correlated, samples may overlap, sample and candidate selection differ; empirical frequency is not a calibrated p-value.'}
for threshold in [0,2,3,4]:
 v=[x for x in records.values() if abs(x[1])>=threshold];zs=[x[0] for x in v];result['groups'][str(threshold)]={'n':len(v),'median_abs_z_difference':statistics.median(abs(x) for x in zs),'mean_squared_z_difference':sum(x*x for x in zs)/len(zs),'fraction_abs_z_gt_1_96':sum(abs(x)>1.95996398454 for x in zs)/len(zs),'fraction_opposite_point_signs':sum(u*t<0 for z,zm,u,t in v)/len(v)}
for c in ca:
 z=float(c['TwinGene_minus_UKB'])/float(c['difference_se_assuming_independence']);meta_z=abs(float(c['stored_meta_beta'])/math.sqrt(1/(1/float(c['UKB_se'])**2+1/float(c['TwinGene_se'])**2)))
 d={'z_difference':z,'meta_z_approx_unrounded_se':meta_z,'empirical_fraction_abs_z_at_least_candidate':{}}
 for threshold in [0,2,3,4]:
  v=[x for x in records.values() if abs(x[1])>=threshold];d['empirical_fraction_abs_z_at_least_candidate'][str(threshold)]=sum(abs(x[0])>=abs(z) for x in v)/len(v)
 result['candidates'][c['variant']]=d
p=pathlib.Path('outputs/fi_cohort_heterogeneity_reference_context.json');p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
