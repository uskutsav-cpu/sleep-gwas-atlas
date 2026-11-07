import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;meta=json.loads((p/'query_receipt.json').read_text())
with (p/'selected_variants.tsv').open(newline='') as h:rows=list(csv.DictReader(h,delimiter='\t'))
assert len(rows)==3 and [r['variant'] for r in rows]==['rs62264764','rs77960','rs1476535']; maxp=0
for row,q in zip(rows,meta['query_records']):
 assert row['variant']==q['variant'] and q['matches']==1 and row['status']=='MATCHED'
 matches=[x.split('\t') for x in q['raw_rows'] if len(x.split('\t'))==12 and x.split('\t')[8]==row['variant']]
 assert len(matches)==1; a=matches[0]; beta=float(a[4])*(1 if a[2]==row['effect_allele'] else -1)
 assert math.isclose(beta,float(row['symptom_beta_aligned']),rel_tol=1e-12)
 assert math.isclose(float(a[7]),float(row['symptom_p']),rel_tol=1e-12)
 assert math.isclose(float(row['symptom_p_times3']),min(1,3*float(a[7])),rel_tol=1e-12)
 wald=math.erfc(abs(float(a[4])/float(a[5]))/math.sqrt(2));maxp=max(maxp,abs(wald-float(a[7])))
status={'rows_verified':len(rows),'max_absolute_wald_vs_source_p_difference':maxp,'selected_variants_sha256':hashlib.sha256((p/'selected_variants.tsv').read_bytes()).hexdigest(),'note':'Orientation and arithmetic checks only; no full-file checksum, independent cohort audit, or validation of case phenotype.'};(p/'verification.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status,indent=2))
