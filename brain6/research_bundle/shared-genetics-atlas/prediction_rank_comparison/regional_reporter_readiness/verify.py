import collections,csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;rows=list(csv.DictReader((p/'regional_rows.tsv').open(),delimiter='\t'));rec=json.loads((p/'receipt.json').read_text());assert rec['zip_crc_verified'] and hashlib.sha256((p/'regional_rows.tsv').read_bytes()).hexdigest()==rec['regional_sha256']
for r in rows:
 s=float(r['Log2Skew']);se=float(r['Skew_SE']);f=10**(-float(r['Skew_logFDR']));pa=10**(-float(r['A_logPadj_BF']));pb=10**(-float(r['B_logPadj_BF']))
 active=(float(r['A_log2FC'])>0 and pa<.05) or (float(r['B_log2FC'])>0 and pb<.05)
 assert (r['proxy_active']=='True')==active and (r['proxy_hit']=='True')==(active and f<.05) and (r['proxy_near_zero']=='True')==(active and abs(s)+1.96*se<.25)
for reg in ('chr3','chr5'):
 a=[r for r in rows if r['region']==reg];q=rec['summary'][reg];assert len(a)==q['rows'] and sum(r['proxy_hit']=='True' for r in a)==q['proxy_hit'] and sum(r['proxy_near_zero']=='True' for r in a)==q['proxy_near_zero']
(p/'verification.json').write_text(json.dumps(dict(status='PASS',checks='input hash, independent p-value transformations and proxy counts',regional_sha256=rec['regional_sha256']),indent=2)+'\n');print('PASS: proxy thresholds and group counts independently recalculated')
