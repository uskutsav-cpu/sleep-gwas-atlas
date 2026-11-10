import csv,math,json,hashlib,pathlib
p=pathlib.Path('outputs/shared22_eight_source_evidence_matrix.csv');rows=list(csv.DictReader(p.open())); assert len(rows)==22 and len({r['variant'] for r in rows})==22
names=['mak_FI','mak_FP','hfrs','finngen_insomnia'];out={}
for n in names:
 diffs=[]; cierrors=[]
 for row in rows:
  b=float(row[n+'_beta_aligned']);s=float(row[n+'_se']);reported=float(row[n+'_p']);assert math.isfinite(b) and math.isfinite(s) and s>0 and 0<reported<=1
  p_wald=math.erfc(abs(b/s)/math.sqrt(2));assert p_wald>0
  diffs.append((abs(math.log10(p_wald/reported)),row['variant']))
  lo,hi=float(row[n+'_ci95_low']),float(row[n+'_ci95_high']);crit=1.959963984540054
  cierrors.extend([abs(lo-(b-crit*s)),abs(hi-(b+crit*s))])
 maximum,variant=max(diffs);out[n]={'n':len(rows),'max_absolute_log10_p_ratio':maximum,'maximum_relative_p_difference':max(abs(10**(math.log10(math.erfc(abs(float(r[n+'_beta_aligned'])/float(r[n+'_se']))/math.sqrt(2))/float(r[n+'_p'])))-1) for r in rows),'max_diff_variant':variant,'max_ci_absolute_error':max(cierrors)}
 assert maximum<.001 and max(cierrors)<1e-11
j={'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'results':out,'scope':'Arithmetic check of stored aligned betas, SE, p and Wald intervals only. Approximate Wald p-values can differ from source tests; no source identity, causal, independence, or validity conclusions.'}
pathlib.Path('outputs/shared22_numerical_consistency_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
