#!/usr/bin/env python3
"""Read-only structural forensic pass; never repairs or admits corrupted data."""
import collections,csv,gzip,json,math,pathlib
root=pathlib.Path('outputs/shared-genetics-atlas/placo_diagnostic');base=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support')
result=root/'FFS_screen.tsv.gz';counts=collections.Counter();first={};seen=0;err=''
def flag(kind,i,detail):
 counts[kind]+=1
 if kind not in first:first[kind]={'row_1based':i,'detail':str(detail)[:180]}
try:
 with gzip.open(base/'S.sumstats.gz','rt') as sf,gzip.open(base/'FFS.sumstats.gz','rt') as tf,gzip.open(result,'rt') as rf:
  s=csv.DictReader(sf,delimiter='\t');t=csv.DictReader(tf,delimiter='\t');r=csv.DictReader(rf,delimiter='\t')
  for i,(a,b,c) in enumerate(zip(s,t,r,strict=True),1):
   seen=i
   if None in c or any(v is None for v in c.values()):flag('malformed_columns',i,c)
   if tuple(a.get(k) for k in ('SNP','A1','A2'))!=tuple(c.get(k) for k in ('SNP','A1','A2')):flag('input_key_mismatch',i,c.get('SNP'))
   if c.get('status') not in ('TESTED','NUMERICAL_FAILURE','EXCLUDED_EXTREME_Z'):flag('invalid_status',i,c.get('status'))
   try:
    za,zb=float(a['Z']),float(b['Z'])
    if float(c['Z_insomnia'])!=za or float(c['Z_frailty'])!=zb:flag('input_z_mismatch',i,c.get('SNP'))
    for z,col in [(za,'P_insomnia'),(zb,'P_frailty')]:
     if not math.isclose(math.erfc(abs(z)/math.sqrt(2)),float(c[col]),rel_tol=2e-12,abs_tol=1e-300):flag('marginal_p_mismatch',i,(c.get('SNP'),col))
    if c.get('status')=='TESTED' and not 0<float(c['P_PLACO_plus'])<=1:flag('invalid_placo_p',i,c.get('SNP'))
   except (KeyError,TypeError,ValueError):flag('unparseable_numeric',i,c.get('SNP'))
except Exception as e:err=type(e).__name__+': '+str(e)
report={'rows_streamed_before_read_error':seen,'structural_mismatch_counts':dict(counts),'first_mismatch':first,'terminal_read_error':err,'interpretation':'Structural checks do not validate PLACO p-values even where fields appear plausible. CRC failure invalidates the complete file; retain it only as forensic evidence and use the full rerun.'}
(root/'ffs_recovery'/'corruption_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
