import gzip,csv,math,collections,hashlib,json,pathlib
p=pathlib.Path(__file__).resolve().parent
root=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support')
expected={'S':'7c0049735306556f46ab3fae3a217a33b172cd25841855061bd9e5d8c95f8bf4','FI':'4d2409b90e557ce12e05e5002822989e5be2db916ceafc727194ecfdab705e4e','FFS':'a16ffd1fb2073dd969416897aca3a81af6dafc57e3b9436a6aa208df901902e6'}
results=[]
for name in ('S','FI','FFS'):
 path=root/(name+'.sumstats.gz');h=hashlib.sha256(path.read_bytes()).hexdigest();assert h==expected[name]
 counts=collections.Counter();n=0;thresholds=(.05,.01,5e-8);uncertain={v:0 for v in thresholds};extreme={v:0 for v in (3,5,6)}
 with gzip.open(path,'rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   z=float(x['Z']);assert math.isfinite(z);v=round(abs(z),3);counts[v]+=1;n+=1
   for cutoff in thresholds:
    upper=math.erfc(max(0,v-.0005)/math.sqrt(2));lower=math.erfc((v+.0005)/math.sqrt(2))
    if lower<=cutoff<upper:uncertain[cutoff]+=1
   for zcut in extreme:
    if v>=zcut:extreme[zcut]+=1
 assert n==1138262
 tied=sum(k for k in counts.values() if k>1)
 results.append({'trait':name,'rows':n,'unique_abs_z':len(counts),'rows_in_ties':tied,'fraction_in_ties':tied/n,'max_abs_z':max(counts),'z_at_least_3':extreme[3],'z_at_least_5':extreme[5],'z_at_least_6':extreme[6],'threshold_ambiguous_0_05':uncertain[.05],'threshold_ambiguous_0_01':uncertain[.01],'threshold_ambiguous_5e_8':uncertain[5e-8],'sha256':h})
with (p/'precision_summary.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
(p/'receipt.json').write_text(json.dumps({'source_rows':1138262,'source_files_sha256':expected,'assumption':'Z rounded to three decimals, independent rounding interval +/- 0.0005; not an error bound for upstream p-to-Z or association model differences','script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print(results)
