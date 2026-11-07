import csv,json,hashlib,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1]
def read(q):
 with q.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(p/'comparison.tsv');b={x['variant']:x for x in read(root/'replication/childhood_symptoms/selected_variants.tsv')};assert len(a)==3 and len(b)==3
err=0
for r in a:
 assert r['status']=='MATCHED';s=float(r['sleep_p']);q=float(b[r['variant']]['symptom_p']);assert float(r['symptom_p'])==q;assert float(r['joint_max_p'])==max(s,q);assert float(r['joint_times3'])==min(1,3*max(s,q));w=math.erfc(abs(float(r['sleep_beta_aligned'])/float(r['sleep_se']))/math.sqrt(2));err=max(err,abs(w-s))
res={'rows_checked':3,'joint_max_and_multiplication':'passed','maximum_sleep_wald_vs_source_p_abs_difference':err,'comparison_sha256':hashlib.sha256((p/'comparison.tsv').read_bytes()).hexdigest(),'note':'Arithmetic consistency only; not correction for selection, causal variant or sample-independence validation.'};(p/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
