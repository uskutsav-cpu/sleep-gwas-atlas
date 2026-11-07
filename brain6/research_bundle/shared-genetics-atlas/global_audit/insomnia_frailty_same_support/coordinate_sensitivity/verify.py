import csv,hashlib,json,pathlib,re
p=pathlib.Path(__file__).resolve().parent;receipt=json.loads((p/'support_receipt.json').read_text());assert receipt['excluded']==537 and receipt['remaining']==1134416
with (p.parent/'correlations.tsv').open() as f:old={x['frailty_trait']:x for x in csv.DictReader(f,delimiter='\t')}
rows=[]
for t in ('FI','FFS','HFRS'):
 text=(p/f'{t}_insomnia.log').read_text();assert '1134416 SNPs with valid alleles.' in text and 'Analysis finished at' in text
 m=re.search(r'Genetic Correlation:\s*([0-9.]+)\s*\(([0-9.]+)\)\s*Z-score:\s*([0-9.]+)\s*P:\s*([0-9.eE+-]+)',text);assert m
 rg,se,z,prob=map(float,m.groups());rows.append({'trait':t,'full_rg':old[t]['rg'],'full_se':old[t]['se'],'restricted_rg':rg,'restricted_se':se,'restricted_p':prob,'difference_rounded':round(rg-float(old[t]['rg']),4),'snps':1134416})
with (p/'results.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(x):
 h=hashlib.sha256()
 with x.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
for name,expected in receipt['outputs'].items():assert sha(pathlib.Path(name))==expected
(p/'verification.json').write_text(json.dumps({'status':'PASS','all_three_completed_logs':True,'input_hashes_match_support_receipt':True,'logs_sha256':{t:sha(p/f'{t}_insomnia.log') for t in ('FI','FFS','HFRS')},'output_sha256':sha(p/'results.tsv'),'limitation':'Support exclusion and results are post hoc and conditional on the chain; does not establish validity of excluded markers.'},indent=2)+'\n');print('verified completed outputs and all four restricted input hashes')
