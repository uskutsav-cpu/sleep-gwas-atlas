#!/usr/bin/env python3
"""Cross-check coverage of HFRS threshold crossings using a different filtered-first join."""
import csv,gzip,collections,json,math,pathlib
p=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation'); o=pathlib.Path('outputs/shared-genetics-atlas/four_trait_marginal');cut=.05/(4*1138262)
source=p/'filtered/HFRS_HM3.tsv'; support=p/'common_support/S.sumstats.gz'; eligible_ids=set()
with source.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  try:v=float(r['P'])
  except ValueError:continue
  if math.isfinite(v) and 0<v<=cut:eligible_ids.add(r['SNP'])
records=collections.defaultdict(list)
with source.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['SNP'] in eligible_ids:records[r['SNP']].append(r)
observed={};present=set();n=0
with gzip.open(support,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1; rows=records.get(r['SNP'],[])
  if r['SNP'] in eligible_ids:present.add(r['SNP'])
  if len(rows)!=1:continue
  x=rows[0];a=(r['A1'],r['A2']);b=(x['A1'],x['A2'])
  if len(set(a))!=2 or any(len(z)!=1 or z not in 'ACGT' for z in a+b) or set(a)!=set(b):continue
  try:pv=float(x['P']);beta=float(x['BETA'])
  except ValueError:continue
  if not (math.isfinite(pv) and 0<pv<=cut and math.isfinite(beta)):continue
  observed[r['SNP']]=pv
assert n==1138262
with (o/'crossings.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
reported={r['SNP']:float(r['HFRS_p']) for r in rows if r['HFRS_status']=='MATCHED' and 'HFRS' in r['threshold_crossing'].split('+')}
assert observed==reported, (len(observed),len(reported),set(observed)^set(reported))
check={'HFRS_source_threshold_snp_ids':len(eligible_ids),'HFRS_threshold_eligible_on_common_support':len(observed),'matches_main_screen_exactly':True,'threshold_ids_absent_from_common_support':sorted(eligible_ids-present),'threshold_ids_present_but_ineligible':sorted(present-set(observed)),'scope':'Independent threshold-first join verifies no HFRS threshold candidate was missed on the support, under the same allele rules; source, chromosome/build and association model validity are not independently established.'}
(o/'independent_hfrs_check.json').write_text(json.dumps(check,indent=2)+'\n');print(json.dumps(check,indent=2))
