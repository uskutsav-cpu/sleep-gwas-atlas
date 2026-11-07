import collections,csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;receipt=json.loads((p/'GCST90624046_receipt.json').read_text());assert receipt['md5']==receipt['publisher_md5'] and receipt['bytes']==420712217
with (p/'GCST90624046_crosswalk.tsv').open() as f:a={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
with (p/'indexed/GCST90624046_extended.tsv').open() as f:b={x['selected_rsid']:x for x in csv.DictReader(f,delimiter='\t')}
assert len(a)==len(b)==87 and a.keys()==b.keys();out=[]
for rs,x in a.items():
 y=b[rs];assert y['lift_status']=='UNIQUE_+';assert x['status']==y['status'],(rs,x['status'],y['status'])
 row={'rsid':rs,'original_status':x['status'],'harmonized_status':y['status'],'effect_allele_flip':'NA','value_comparison':'NOT_APPLICABLE'}
 if x['status']=='UNIQUE_MATCH':
  assert {x['effect_allele'],x['other_allele']}=={y['effect_allele'],y['other_allele']};flip=1 if x['effect_allele']==y['effect_allele'] else -1
  assert math.isclose(float(x['beta'])*flip,float(y['beta']),rel_tol=1e-10,abs_tol=1e-14)
  f=float(x['effect_allele_frequency']);f=f if flip==1 else 1-f;assert math.isclose(f,float(y['effect_allele_frequency']),rel_tol=1e-10,abs_tol=1e-14)
  for k in ('standard_error','p_value','Q_metric','Q_metric_pvalue'):assert math.isclose(float(x[k]),float(y[k]),rel_tol=1e-10,abs_tol=1e-14),(rs,k)
  row.update(effect_allele_flip=flip,value_comparison='PASS')
 out.append(row)
with (p/'reconciliation.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary={'selected':len(out),'status_counts':dict(collections.Counter(r['original_status'] for r in out)),'matched_value_comparisons_passed':sum(r['value_comparison']=='PASS' for r in out),'matched_allele_flips':dict(collections.Counter(str(r['effect_allele_flip']) for r in out if r['value_comparison']=='PASS')),'original_release_md5_verified':True,'harmonized_full_file_md5_verified':False,'limitation':'Common absences show absence from two filtered releases, not the unfiltered GWAS or no association; same-study consistency is not replication.'};(p/'reconciliation_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
