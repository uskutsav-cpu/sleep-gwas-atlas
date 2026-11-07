import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1]
def read(q):
 with q.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
x=read(p/'six_leads.tsv');cand=[r for r in read(root/'evidence/BRAIN6_FINAL_CANDIDATES.tsv') if r['pair_id']=='insomnia__adhd'];rec=json.loads((p/'receipt.json').read_text());assert len(x)==len(cand)==6 and {r['variant'] for r in x}=={r['lead_variants'] for r in cand}
errors=[]
for r in x:
 if r['status']!='MATCHED':continue
 q=[z for z in rec['query_records'] if z['variant']==r['variant']];assert len(q)==1 and q[0]['sleep_matches']==q[0]['symptom_matches']==1
 vals=[z.split('\t') for z in q[0]['symptom_raw_rows'] if len(z.split('\t'))==12 and z.split('\t')[8]==r['variant']];assert len(vals)==1;v=vals[0]
 aligned=float(v[4])*(1 if v[2]==r['effect_allele'] else -1);assert math.isclose(aligned,float(r['symptom_beta_aligned']),rel_tol=1e-12)
 a,b=float(r['sleep_p']),float(r['symptom_p']);assert math.isclose(float(r['max_p']),max(a,b),rel_tol=1e-12) and math.isclose(float(r['max_p_times_six']),min(1,6*max(a,b)),rel_tol=1e-12)
 w=math.erfc(abs(float(r['sleep_beta_aligned'])/float(r['sleep_se']))/math.sqrt(2));errors.append(abs(w-a))
out={'leads':len(x),'matched':sum(r['status']=='MATCHED' for r in x),'not_estimable':sum(r['status']!='MATCHED' for r in x),'maximum_sleep_wald_vs_source_p_difference':max(errors),'file_sha256':hashlib.sha256((p/'six_leads.tsv').read_bytes()).hexdigest(),'note':'Arithmetic, frozen-family membership and symptom alignment only; no scientific calibration or full-family selection control.'};(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
