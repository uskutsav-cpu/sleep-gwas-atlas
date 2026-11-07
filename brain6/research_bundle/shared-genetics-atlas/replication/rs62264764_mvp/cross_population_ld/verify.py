import csv,json,math,pathlib,hashlib
p=pathlib.Path(__file__).resolve().parent; rows=list(csv.DictReader((p/'group_ld.tsv').open(),delimiter='\t'));rec=json.loads((p/'receipt.json').read_text())
assert len(rows)==5 and {r['group'] for r in rows}=={'AFR','AMR','EAS','EUR','SAS'}
for r in rows:
 n=int(r['haplotypes']);a,b,c,d=(int(r[k]) for k in ('n00','n01','n10','n11'));assert a+b+c+d==n==2*int(r['complete_n'])
 f1=(c+d)/n;f2=(b+d)/n;v=((a*d-b*c)**2)/((a+b)*(c+d)*(a+c)*(b+d))
 assert math.isclose(f1,float(r['freq1']),abs_tol=1e-14) and math.isclose(f2,float(r['freq2']),abs_tol=1e-14) and math.isclose(v,float(r['r2']),abs_tol=1e-12)
 assert int(r['bootstrap_valid'])==2000
assert next(r for r in rows if r['group']=='EUR')['n11']=='165'
(p/'verification.json').write_text(json.dumps(dict(status='PASS',method='independent haplotype-count cross-product formula, frequencies and group coverage; no bootstrap regeneration',input_sha256=hashlib.sha256((p/'group_ld.tsv').read_bytes()).hexdigest()),indent=2)+'\n')
print('PASS: independently recomputed five frequencies and r² values from haplotype counts')
