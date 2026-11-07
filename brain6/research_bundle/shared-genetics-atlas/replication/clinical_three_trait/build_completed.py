import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def rows(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
old=p/'results.tsv'; completion=p/'rs2582894_completion.tsv'; receipt=p/'completion_receipt.json'; prior=rows(old); extra=rows(completion); cert=json.loads(receipt.read_text()); assert len(prior)==7 and len(extra)==1 and cert['source_md5']=='38d9434f50a4493def9d2c9d1cc1b4b3' and cert['source_rows']==21294560
assert cert['input_hashes'][str(old)]==sha(old)
for path,hashval in cert['input_hashes'].items():assert sha(Path(path))==hashval
x=extra[0]; assert x['variant']=='rs2582894'; z=[r for r in prior if r['variant']==x['variant']];assert len(z)==1 and z[0]['status']=='not_estimable';z=z[0]
for k in ('grch38','mvp_insomnia_p','finngen_adhd_p','finngen_hfrs_p','adhd_beta_alt','hfrs_beta_alt','max_p','max_p_times_7'):assert str(cert['output'][k])==x[k]
assert str(cert['output']['mvp_beta_alt'])==x['mvp_beta_alt']
z.update({'candidate_status':'completed_from_full_hfrs','grch38':x['grch38'],'mvp_insomnia_p':x['mvp_insomnia_p'],'finngen_adhd_p':x['finngen_adhd_p'],'finngen_hfrs_p':x['finngen_hfrs_p'],'adhd_beta_alt':x['adhd_beta_alt'],'insomnia_beta_alt':x['mvp_beta_alt'],'hfrs_beta_alt':x['hfrs_beta_alt'],'max_p':x['max_p'],'max_p_times_7':x['max_p_times_7'],'status':'matched'})
for row in prior:
 if row['status']=='matched':
  ps=[float(row[k]) for k in ('mvp_insomnia_p','finngen_adhd_p','finngen_hfrs_p')];assert all(math.isfinite(q) and 0<=q<=1 for q in ps);m=max(ps); assert math.isclose(float(row['max_p']),m,abs_tol=1e-15) and math.isclose(float(row['max_p_times_7']),min(1,7*m),abs_tol=1e-15)
assert sum(x['status']=='matched' for x in prior)==6 and len({x['variant'] for x in prior})==7
with (p/'completed_results.tsv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=prior[0].keys(),delimiter='\t');w.writeheader();w.writerows(prior)
(p/'completed_table_receipt.json').write_text(json.dumps({'input_sha256':{str(f):sha(f) for f in (old,completion,receipt)},'output_sha256':sha(p/'completed_results.tsv'),'slots':len(prior),'complete':6,'missing':1,'warning':'arithmetic and provenance checks do not validate source data or establish independent replication'},indent=2)+'\n')
print('Rebuilt and checked six complete three-trait rows, one missing row')
