import csv,gzip,hashlib,json,math
from pathlib import Path

here=Path(__file__).resolve().parent
out=Path('/Volumes/Brain6Work/physical-function-source-sensitivity-20261028')
out.mkdir(exist_ok=False)
paths={'walking':Path('/Volumes/Brain6Work/neale_924_walking.ldsc.bgz'),
       'grip':Path('/Volumes/Brain6Work/neale_46irnt_grip.ldsc.bgz'),
       'registry':Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/S.sumstats.gz'),
       'questionnaire':Path('/Volumes/Brain6Work/brain6-frailty-paired/insomnia.sumstats.gz')}
files={k:gzip.open(v,'rt') for k,v in paths.items()}
readers={k:csv.DictReader(f,delimiter='\t') for k,f in files.items()}
outputs={k:gzip.open(out/f'{k}.sumstats.gz','wt') for k in paths}
writers={k:csv.DictWriter(f,fieldnames=['SNP','A1','A2','N','Z'],delimiter='\t',lineterminator='\n') for k,f in outputs.items()}
for w in writers.values():w.writeheader()
last={k:'' for k in paths}; counts={k:0 for k in paths}; match=0; skipped=0
def advance(k):
    row=next(readers[k],None)
    if row is not None:
        assert row['SNP']>last[k],(k,row['SNP'],last[k])
        last[k]=row['SNP'];counts[k]+=1
        assert math.isfinite(float(row['Z'])) and 0<float(row['N'])<1e9
    return row
rows={k:advance(k) for k in paths};comp=str.maketrans('ACGT','TGCA')
while all(x is not None for x in rows.values()):
    hi=max(r['SNP'] for r in rows.values())
    for k in paths:
        while rows[k] is not None and rows[k]['SNP']<hi:rows[k]=advance(k)
    if not all(x is not None for x in rows.values()):break
    if len({x['SNP'] for x in rows.values()})!=1:continue
    alleles=[(r['A1'],r['A2']) for r in rows.values()]; a,b=alleles[0]
    good=len(a)==len(b)==1 and a in 'ACGT' and b in 'ACGT' and a!=b and {a,b} not in ({'A','T'},{'C','G'})
    good=good and all({c,d}=={a,b} or {c.translate(comp),d.translate(comp)}=={a,b} for c,d in alleles)
    if good:
        for k,w in writers.items():w.writerow({x:rows[k][x] for x in w.fieldnames})
        match+=1
    else:skipped+=1
    rows={k:advance(k) for k in paths}
for k in paths:
    while rows[k] is not None:rows[k]=advance(k)
    files[k].close();outputs[k].close()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
receipt={'source_paths':{k:str(v) for k,v in paths.items()},'source_sha256':{k:sha(v) for k,v in paths.items()},
         'rows_read':counts,'all_source_streams_finished':True,'common_valid_nonpalindromic':match,
         'shared_but_invalid_alleles':skipped,'output_paths':{k:str(out/f'{k}.sumstats.gz') for k in paths},
         'output_sha256':{k:sha(out/f'{k}.sumstats.gz') for k in paths}}
(here/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
