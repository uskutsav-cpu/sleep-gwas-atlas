import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1]
def read(q):
 with q.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
x=read(p/'six_leads.tsv');r=json.loads((p/'receipt.json').read_text());prior=json.loads((root/'replication/six_lead_joint/receipt.json').read_text());assert len(x)==6 and r['mvp_rows']==19703815
for q in x:
 a=[z for z in r['audit'] if z['variant']==q['variant']];b=[z for z in prior['query_records'] if z['variant']==q['variant']];assert len(a)==len(b)==1 and a[0]['mvp_exact_matches']==a[0]['symptom_exact_matches']==1
 mv=a[0]['mvp_exact_rows'][0];s=[z.split('\t') for z in b[0]['symptom_raw_rows'] if len(z.split('\t'))==12 and z.split('\t')[8]==q['variant'] and q['effect_allele'] in z.split('\t')[2:4]];assert len(s)==1
 logodds=math.log(float(mv['odds_ratio']))*(1 if mv['effect_allele']==q['effect_allele'] else -1);beta=float(s[0][4])*(1 if s[0][2]==q['effect_allele'] else -1)
 assert math.isclose(logodds,float(q['mvp_logor_aligned']),abs_tol=1e-14) and math.isclose(beta,float(q['symptom_beta_aligned']),abs_tol=1e-14)
 assert float(q['mvp_p'])==float(mv['p_value']) and float(q['symptom_p'])==float(s[0][7]);assert float(q['joint_max_p'])==max(float(q['mvp_p']),float(q['symptom_p'])) and math.isclose(float(q['joint_times_six']),min(1,6*float(q['joint_max_p'])),abs_tol=1e-14)
out={'checked':len(x),'all_exact_source_matches':True,'max_p_and_allele_alignment':'passed','sha256':hashlib.sha256((p/'six_leads.tsv').read_bytes()).hexdigest(),'note':'Independent arithmetic check, not scientific replication, sample independence or source-model validation.'};(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
