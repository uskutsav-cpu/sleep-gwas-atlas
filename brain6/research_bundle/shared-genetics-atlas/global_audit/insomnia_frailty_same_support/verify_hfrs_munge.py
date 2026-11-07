import csv,gzip,pathlib,hashlib,json,math,collections
import numpy as np
from scipy.special import ndtri
p=pathlib.Path(__file__).resolve().parent;src=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered/HFRS_HM3.tsv');out=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/HFRS.sumstats.gz');log=src.parent.parent/'munged/HFRS.log';comp=str.maketrans('ACGT','TGCA')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
records={};duplicates=set();counts=collections.Counter()
with src.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  counts['source_rows']+=1;rs=r['SNP']
  if rs in records:duplicates.add(rs);continue
  records[rs]=(r['A1'],r['A2'],r['BETA'],r['P'],r['N'])
probs=[];signs=[];zs=[];bad=[]
with gzip.open(out,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  counts['selected_rows']+=1;s=records.get(r['SNP'])
  if not s or r['SNP'] in duplicates:counts['missing_or_duplicate']+=1;bad.append((r['SNP'],'missing_or_duplicate'));continue
  a,b,beta,pvalue,n=s;x,y=r['A1'],r['A2']
  if any(len(v)!=1 or v not in 'ACGT' for v in (a,b,x,y)) or ''.join(sorted(a+b)) in ('AT','CG'):
   counts['invalid_alleles']+=1;bad.append((r['SNP'],'invalid_alleles'));continue
  if (a,b)==(x,y):flip=1;mode='direct'
  elif (a.translate(comp),b.translate(comp))==(x,y):flip=1;mode='complement'
  elif (b,a)==(x,y):flip=-1;mode='swap'
  elif (b.translate(comp),a.translate(comp))==(x,y):flip=-1;mode='complement_swap'
  else:counts['allele_mismatch']+=1;bad.append((r['SNP'],'allele_mismatch'));continue
  try:
   beta=float(beta);prob=float(pvalue);N=float(n);Z=float(r['Z']);assert all(math.isfinite(v) for v in (beta,prob,N,Z)) and 0<prob<=1 and N>0 and float(r['N'])==N and (beta!=0 or prob==1)
  except (ValueError,AssertionError):counts['invalid_numeric_or_n']+=1;bad.append((r['SNP'],'invalid_numeric_or_n'));continue
  probs.append(prob);signs.append((1 if beta>=0 else -1)*flip);zs.append(Z);counts[mode]+=1
prob=np.asarray(probs);expected=-ndtri(prob/2)*np.asarray(signs);observed=np.asarray(zs);diff=np.abs(expected-observed);invalid=~np.isfinite(expected);mismatch=(diff>0.000501)|invalid
counts['signed_z_mismatch']=int(mismatch.sum());counts['z_compared']=len(diff);assert counts['selected_rows']==1134953 and not bad and not mismatch.any(),(counts,bad[:10],float(np.max(diff)))
r={'source_path':str(src),'source_sha256':sha(src),'final_hfrs_path':str(out),'final_hfrs_sha256':sha(out),'munging_log_path':str(log),'munging_log_sha256':sha(log),'counts':dict(counts),'max_absolute_z_rounding_difference':float(np.max(diff)),'status':'PASS','limitations':'Filtered-source-to-final numerical/allele transformation only; source effect allele, cohort overlap, sample size and ancestry are not independently validated.'};(p/'hfrs_munging_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['counts'],indent=2));print('max z diff',r['max_absolute_z_rounding_difference'])
