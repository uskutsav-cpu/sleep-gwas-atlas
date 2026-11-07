import csv, gzip, hashlib, json, math
from pathlib import Path

p = Path(__file__).resolve().parent
out = Path('/Volumes/Extreme SSD/brain6-physical-factor-20261101')
paths = {x: out / (x + '.bgz') for x in ('walking', 'grip', 'activity')}
paths['insomnia'] = Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz')
expected = {'walking': '0dcc63144c15b10dba6719314fd294439ce11293c13010cdec71f4e9b274535f',
            'grip': '6054617706c15b99be4df40ef731a1e077dcb729806d53bd403235935338b987',
            'activity': 'fffebe514901e7c09b242365f3a327d0a73d4df6729d64fb903f45a4eb7ad5be',
            'insomnia': '3d7073182abbddfbd83eefbd1187ac999a5d49321ec982d19c851ac1027826c6'}
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while b := f.read(1048576): h.update(b)
    return h.hexdigest()
observed = {k: sha(v) for k, v in paths.items()}
if observed != expected: raise ValueError(('source hash mismatch', observed))
files = {k: gzip.open(v, 'rt') for k,v in paths.items()}
readers = {k: csv.DictReader(v, delimiter='\t') for k,v in files.items()}
for k,r in readers.items():
    if set(r.fieldnames) != {'SNP','A1','A2','N','Z'}: raise ValueError((k,r.fieldnames))
outputs = {k: gzip.open(out / (k + '.common.sumstats.gz'), 'xt') for k in paths}
writers = {k: csv.DictWriter(v, fieldnames=['SNP','A1','A2','N','Z'], delimiter='\t', lineterminator='\n') for k,v in outputs.items()}
for w in writers.values(): w.writeheader()
counts = {k:0 for k in paths}; last={k:'' for k in paths}; common=bad=0
def advance(k):
    row=next(readers[k],None)
    if row is not None:
        if row['SNP'] <= last[k]: raise ValueError(('unsorted or duplicate', k, row['SNP']))
        last[k]=row['SNP']; counts[k]+=1
        if not math.isfinite(float(row['Z'])) or not 0<float(row['N'])<1e9: raise ValueError(('invalid', k, row))
    return row
rows={k:advance(k) for k in paths}; complement=str.maketrans('ACGT','TGCA')
while all(x is not None for x in rows.values()):
    highest=max(x['SNP'] for x in rows.values())
    for k in paths:
        while rows[k] is not None and rows[k]['SNP'] < highest: rows[k]=advance(k)
    if not all(x is not None for x in rows.values()): break
    if len({x['SNP'] for x in rows.values()}) != 1: continue
    a,b=rows['walking']['A1'],rows['walking']['A2']
    good=len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and {a,b} not in ({'A','T'},{'C','G'})
    good=good and all({r['A1'],r['A2']} == {a,b} or {r['A1'].translate(complement),r['A2'].translate(complement)} == {a,b} for r in rows.values())
    if good:
        for k,w in writers.items(): w.writerow({x: rows[k][x] for x in w.fieldnames})
        common+=1
    else: bad+=1
    rows={k:advance(k) for k in paths}
for k in paths:
    while rows[k] is not None: rows[k]=advance(k)
    files[k].close(); outputs[k].close()
receipt={'source_sha256':observed,'input_rows':counts,'common_nonpalindromic_rows':common,'excluded_shared':bad,
         'output_sha256':{k:sha(out/(k+'.common.sumstats.gz')) for k in paths},'output_directory':str(out)}
(p/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps(receipt,indent=2))
