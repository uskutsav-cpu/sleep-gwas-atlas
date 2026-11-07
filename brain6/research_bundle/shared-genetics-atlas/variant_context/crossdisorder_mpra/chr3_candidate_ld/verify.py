import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;rows=list(csv.DictReader((p/'ld.tsv').open(),delimiter='\t'));rec=json.loads((p/'receipt.json').read_text());assert len(rows)==4 and all(r['status']=='matched' for r in rec['statuses'])
for r in rows:
 a,b,c,d=(int(r[k]) for k in ('n00','n01','n10','n11'));n=a+b+c+d;assert n==2*int(r['n'])
 r2=((a*d-b*c)**2)/((a+b)*(c+d)*(a+c)*(b+d));assert math.isclose(r2,float(r['r2']),abs_tol=1e-12) and int(r['bootstrap_valid'])==2000
 assert math.isclose((c+d)/n,float(r['freq1']),abs_tol=1e-14) and math.isclose((b+d)/n,float(r['freq2']),abs_tol=1e-14)
(p/'verification.json').write_text(json.dumps(dict(status='PASS',check='independent haplotype cross-product r², frequencies, counts and match status; bootstrap not rerun',ld_sha256=hashlib.sha256((p/'ld.tsv').read_bytes()).hexdigest()),indent=2)+'\n')
print('PASS: all four pairwise r² values independently checked')
