import csv,json,math,pathlib,hashlib
p=pathlib.Path(__file__).resolve().parent;rows=list(csv.DictReader((p/'source_rows.tsv').open(),delimiter='\t'));rec=json.loads((p/'receipt.json').read_text());assert len(rows)==4 and len(rec['files'])==4
for row in rows:
 assert row['chr'].removeprefix('chr')=='3' and row['pos']=='117637697' and (row['ref_allele'],row['alt_allele'])==('T','C')
 for name in ('Skew_logP','Skew_logFDR','A_logPadj_BF','B_logPadj_BF'):
  v=10**(-float(row[name]));assert math.isfinite(v) and 0<=v<=1
assert all(x['zip_crc_verified'] and x['matches']==1 for x in rec['files'])
(p/'verification.json').write_text(json.dumps(dict(status='PASS',checks='unique exact position and alleles, four member receipts, transformed p-value ranges',source_rows_sha256=hashlib.sha256((p/'source_rows.tsv').read_bytes()).hexdigest()),indent=2)+'\n')
print('PASS: exact matches, four CRC receipts and p-value transformations')
