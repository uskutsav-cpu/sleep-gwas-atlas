#!/usr/bin/env python3
"""Prepare an allele-oriented comparison AFTER the full stream succeeds.

Does not certify stream completeness. Check shell exit and gzip CRC first.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--verified-gzip',action='store_true',help='Assert that the full source pipeline exited 0 with gzip CRC verified')
a=parser.parse_args()
if not a.verified_gzip: parser.error('Requires verification that the full pipeline exited zero')
def read(name):
 with (p/name).open(newline='') as f:return list(csv.DictReader(f))
raw=read('stage6_bbj_insomnia_candidate_matches.csv')
defs={}
for label,file in [('prototype15','unvalidated_15_variant_score_prototype.csv'),('shared22','unvalidated_22_shared_direction_score.csv')]:
 for x in read(file):
  key=(label,x['variant'])
  if key in defs:raise ValueError(f'duplicate candidate {key}')
  defs[key]=x
by={}
for x in raw:
 key=(x['score'],x['candidate_id'])
 if x['score']=='paper_control':continue
 if key not in defs:raise ValueError(f'unknown candidate {key}')
 if x['exact_pair']=='1':by.setdefault(key,[]).append(x)
if any(len(v)>1 for v in by.values()):raise ValueError('multiple exact allele matches; manual review needed')
reference={x['variant']:x for x in read('stage6_shared22_cross_population_frequencies.csv')}
result=[]
for (score,variant),candidate in defs.items():
 row={'score':score,'variant':variant,'counted_allele':candidate['counted_allele'],'other_allele':candidate['other_allele']}
 matches=by.get((score,variant),[])
 if not matches:
  row.update({'status':'no exact pair in extracted autosomal source'})
 else:
  x=matches[0];counted=row['counted_allele'];other=row['other_allele']
  if {x['Allele1'],x['Allele2']}!={counted,other}:raise ValueError('allele pair mismatch')
  beta=float(x['BETA']); se=float(x['SE']); af=float(x['AF_Allele2']); info=float(x['imputationInfo']); n=float(x['N']); ac=float(x['AC_Allele2'])
  if not math.isfinite(beta) or not (math.isfinite(se) and se>0 and 0<=af<=1 and math.isfinite(info)):raise ValueError('bad numeric field')
  oriented=beta if x['Allele2']==counted else -beta
  freq=af if x['Allele2']==counted else 1-af
  palindrome={counted,other} in ({'A','T'},{'C','G'})
  ref=reference.get(variant) if score=='shared22' else None
  if ref and (ref['counted_allele'],ref['other_allele'])!=(counted,other):raise ValueError('reference allele mismatch '+variant)
  ref_freq=float(ref['EAS_counted_frequency']) if ref else None
  row.update({'reference_EAS_counted_frequency':ref_freq,'absolute_EAS_BBJ_frequency_difference':abs(freq-ref_freq) if ref else None,'status':'exact position and pair','source_id':x['SNPID'],'source_allele1':x['Allele1'],'source_allele2':x['Allele2'],'effect_allele_assumption':'Allele2 per SAIGE output documentation','palindromic':palindrome,'counted_allele_frequency':freq,'imputation_info':info,'allele2_dosage_count':ac,'minor_dosage_count':min(ac,2*n-ac),'minor_dosage_count_under_20':min(ac,2*n-ac)<20,'allele2_frequency_in_cases':x['AF.Cases'],'allele2_frequency_in_controls':x['AF.Controls'],'beta_counted_assuming_saige':oriented,'se':se,'wald_95_low':oriented-1.96*se,'wald_95_high':oriented+1.96*se,'source_p_spa':x['p.value'],'source_p_normal':x['p.value.NA'],'source_n':x['N']})
 result.append(row)
fields=['score','variant','counted_allele','other_allele','status','source_id','source_allele1','source_allele2','effect_allele_assumption','palindromic','counted_allele_frequency','reference_EAS_counted_frequency','absolute_EAS_BBJ_frequency_difference','imputation_info','allele2_dosage_count','minor_dosage_count','minor_dosage_count_under_20','allele2_frequency_in_cases','allele2_frequency_in_controls','beta_counted_assuming_saige','se','wald_95_low','wald_95_high','source_p_spa','source_p_normal','source_n']
out=p/'stage6_bbj_insomnia_oriented_candidates.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(result)
counts={score:{'total':sum(x['score']==score for x in result),'exact':sum(x['score']==score and x['status']=='exact position and pair' for x in result)} for score in ('prototype15','shared22')}
receipt={'candidate_source_sha256':hashlib.sha256((p/'stage6_bbj_insomnia_candidate_matches.csv').read_bytes()).hexdigest(),'counts':counts,'source_allele_interpretation':'conditional on SAIGE Allele2 coding; source-specific strand and multiallelic representations still require review','caveats':'candidate association follow-up only, not score validation; Wald intervals may not correspond to SPA p values'}
(p/'stage6_bbj_insomnia_oriented_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(counts,indent=2))
