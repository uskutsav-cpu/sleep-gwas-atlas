import csv,gzip,hashlib,json,math,re
from pathlib import Path
P=Path('outputs/shared-genetics-atlas/global_audit/brain6_frailty_covariance')
base=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged')
frail=Path('/Volumes/Extreme SSD/sleep-frailty-full-tp-v3-20261009/munged_full_native')
files={'insomnia':base/'insomnia.sumstats.gz','adhd':base/'adhd.sumstats.gz','fi':frail/'frailty.sumstats.gz','ffs':frail/'fried_frailty_score.sumstats.gz'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
old1=json.loads(Path('outputs/shared-genetics-atlas/global_audit/insomnia_adhd_rerun/receipt.json').read_text())['input_frozen_hashes']
old2=json.loads(Path('outputs/shared-genetics-atlas/global_audit/fi_ffs_rerun/receipt.json').read_text())['source_sha256']
for name in ('insomnia','adhd'): assert sha(files[name])==old1[name]
for name in ('fi','ffs'): assert sha(files[name])==old2[str(files[name])]
# Verify the same software and 22 LD-score / M reference files frozen in the earlier Brain6 rerun.
previous=json.loads(Path('outputs/shared-genetics-atlas/global_audit/insomnia_adhd_rerun/receipt.json').read_text())['source_sha256']
reference={k:v for k,v in previous.items() if '/eur_w_ld_chr/' in k or k.endswith('/ldsc/ldsc.py')}
assert len(reference)==45 and all(sha(Path(k))==v for k,v in reference.items())
rows=[]
for name in ['insomnia_frailty','insomnia_fried_frailty_score','adhd_frailty']:
 t=(P/(name+'.log')).read_text(); assert 'Analysis finished at' in t
 n=int(re.search(r'After merging with summary statistics, (\d+) SNPs remain',t).group(1))
 sub=t.split('Genetic Covariance\n')[1]; intercept,se=map(float,re.search(r'Intercept: ([\d.-]+) \(([\d.]+)\)',sub).groups());rg,rgse=map(float,re.search(r'Genetic Correlation: ([\d.-]+) \(([\d.]+)\)',sub).groups());p=float(re.search(r'^P: ([\deE+.-]+)',sub,re.M).group(1))
 rows.append(dict(pair=name,n=n,rg=rg,rg_se=rgse,p=p,bonferroni_3=min(1,3*p),cov_intercept=intercept,cov_intercept_se=se))
with (P/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
# Check local matched FI rows against the native munged signed statistics when an exact nonpalindromic allele match exists.
with open('/Volumes/Brain6Work/chr7-three-trait-input.tsv') as f: local={r['rsid']:r for r in csv.DictReader(f,delimiter='\t')}
matched=aligned=opposite=zero=0; differences=[]; found=set()
with gzip.open(files['fi'],'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  sid=r.get('SNP'); q=local.get(sid)
  if q is None or not r.get('A1') or not r.get('A2') or not r.get('Z'): continue
  if {r['A1'],r['A2']}!={q['effect_allele'],q['other_allele']}:continue
  assert sid not in found;found.add(sid);matched+=1
  z=float(r['Z'])*(1 if r['A1']==q['effect_allele'] else -1); zlocal=float(q['fi_beta'])/float(q['fi_se'])
  differences.append(abs(z-zlocal))
  if z==0 or zlocal==0:zero+=1
  elif z*zlocal>0:aligned+=1
  else:opposite+=1
receipt={'source_sha256':{str(p):sha(p) for p in files.values()},'same_files_as_prior_rechecks':True,'ldsc_and_44_reference_files_match_previous_receipt':True,'results':rows,'local_fi_concordance':{'local_markers':len(local),'matched_valid_alleles':matched,'same_sign':aligned,'opposite_sign':opposite,'zero_sign':zero,'max_abs_difference_native_z_vs_beta_over_se':max(differences) if differences else None,'note':'Descriptive match is not proof of identical underlying source and does not require exact equality between p-derived Z and rounded beta/SE.'},'logs_sha256':{name:sha(P/(name+'.log')) for name in ['insomnia_frailty','insomnia_fried_frailty_score','adhd_frailty']}}
(P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt['local_fi_concordance'],indent=2));print(rows)
