import csv,gzip,hashlib,json
from pathlib import Path

here=Path(__file__).resolve().parent
base=Path('/Volumes/Brain6Work/brain6-frailty-paired')
work=Path('/Volumes/Brain6Work/physical-function-source-sensitivity-20261028')
up=json.loads((here.parent/'prepare_receipt.json').read_text())
prior=json.loads(Path('outputs/shared-genetics-atlas/global_audit/brain6_frailty_covariance/three_frailty_common_support/matching_receipt.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
source={k:base/f'{k}.sumstats.gz' for k in ['FI','FFS','HFRS','insomnia']}
for k in ['FI','FFS','insomnia']: assert sha(source[k])==prior['earlier_input_sha256'][str(source[k])]
assert sha(source['HFRS'])==prior['hfrs_input_sha256']
common=Path(up['output_paths']['questionnaire']);assert sha(common)==up['output_sha256']['questionnaire']
fs={k:gzip.open(v,'rt') for k,v in source.items()};read={k:csv.DictReader(v,delimiter='\t') for k,v in fs.items()}
os={k:gzip.open(work/f'{k}.sumstats.gz','xt') for k in ['FI','FFS','HFRS']}
wr={k:csv.DictWriter(v,fieldnames=['SNP','A1','A2','N','Z'],delimiter='\t',lineterminator='\n') for k,v in os.items()}
for w in wr.values():w.writeheader()
count=selected=0
with gzip.open(common,'rt') as f:
 wanted=csv.DictReader(f,delimiter='\t');target=next(wanted,None)
 for entries in zip(*read.values()):
  rows=dict(zip(read,entries));count+=1
  key=lambda r:(r['SNP'],r['A1'],r['A2'])
  assert len({key(r) for r in rows.values()})==1
  if target is not None and key(rows['insomnia'])==key(target):
   for k,w in wr.items():w.writerow({x:rows[k][x] for x in w.fieldnames})
   selected+=1;target=next(wanted,None)
 assert target is None and selected==up['common_valid_nonpalindromic'] and count==prior['matching_rows']
 assert all(next(it,None) is None for it in read.values())
for f in list(fs.values())+list(os.values()):f.close()
result={'source_sha256':{k:sha(v) for k,v in source.items()},'parent_common_sha256':sha(common),
        'source_rows':count,'selected_rows':selected,'ordered_allele_keys_identical':True,
        'output_sha256':{k:sha(work/f'{k}.sumstats.gz') for k in wr},'work_path':str(work)}
(here/'prepare_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
