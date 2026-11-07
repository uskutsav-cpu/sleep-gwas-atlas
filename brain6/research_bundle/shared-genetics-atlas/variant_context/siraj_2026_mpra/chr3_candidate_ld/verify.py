import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;rows=list(csv.DictReader((p/'ld.tsv').open(),delimiter='\t'));assert {r['rsid'] for r in rows}=={'rs73187947','rs9815484','rs1456193'}
for r in rows:
 a,b,c,d=(int(r[k]) for k in ('n00','n01','n10','n11'));n=a+b+c+d;assert n==2*int(r['n']);v=(a*d-b*c)**2/((a+b)*(c+d)*(a+c)*(b+d));assert math.isclose(v,float(r['r2']),abs_tol=1e-12)
 assert int(r['bootstrap_valid'])==2000
assert abs(float(next(r for r in rows if r['rsid']=='rs1456193')['r2'])-0.821187344472094)<1e-12
q=p/'other_cells';r=list(csv.DictReader((q/'source_rows.tsv').open(),delimiter='\t'));receipt=json.loads((q/'receipt.json').read_text());assert len(r)==4 and len(receipt['files'])==4
for x in r:
 assert x['chr'].removeprefix('chr')=='3' and x['pos']=='117637295' and (x['ref_allele'],x['alt_allele'])==('A','G')
 assert all(0<=10**(-float(x[k]))<=1 for k in ('Skew_logP','Skew_logFDR','A_logPadj_BF','B_logPadj_BF'))
assert all(f['zip_crc_verified'] and f['matches']==1 for f in receipt['files'])
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();(p/'verification.json').write_text(json.dumps(dict(status='PASS',checks='independent LD count formula, prior control, assay allele matching and CRC receipts; bootstrap not regenerated',ld_sha256=sha(p/'ld.tsv'),assay_rows_sha256=sha(q/'source_rows.tsv')),indent=2)+'\n');print('PASS: LD counts, prior control, allele matches and receipt checks')
