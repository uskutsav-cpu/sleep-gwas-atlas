import csv,gzip,hashlib,json
from pathlib import Path

p=Path(__file__).resolve().parent
work=Path('/Volumes/Brain6Work/physical-function-source-sensitivity-20261028')
raw=Path('/Volumes/Brain6Work/ukb_hfrs_public_20261028.sumstats.gz')
stream=json.loads((p/'stream_receipt.json').read_text())
prev=json.loads(Path('outputs/shared-genetics-atlas/component_exclusion/physical_function_ldsc/insomnia_source_sensitivity/prepare_receipt.json').read_text())
fprev=json.loads(Path('outputs/shared-genetics-atlas/component_exclusion/physical_function_ldsc/insomnia_source_sensitivity/frailty_extension/prepare_receipt.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(raw)==stream['output_sha256']
ids=[]
with gzip.open(raw,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):ids.append(r['SNP'])
assert len(ids)==stream['counts']['selected'] and len(ids)==len(set(ids))
out=Path('/Volumes/Brain6Work/ukb-hfrs-comparison-20261028');out.mkdir(exist_ok=False)
sources={'questionnaire':work/'questionnaire.sumstats.gz','registry':work/'registry.sumstats.gz','finngen_hfrs':work/'HFRS.sumstats.gz'}
results={}
for k,source in sources.items():
 expected=fprev['output_sha256']['HFRS'] if k=='finngen_hfrs' else prev['output_sha256'][k]
 assert sha(source)==expected
 match=0; seen=0
 with gzip.open(source,'rt') as f,gzip.open(out/f'{k}.sumstats.gz','xt') as g:
  reader=csv.DictReader(f,delimiter='\t');writer=csv.DictWriter(g,fieldnames=['SNP','A1','A2','N','Z'],delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in reader:
   seen+=1
   if match<len(ids) and row['SNP']==ids[match]:
    writer.writerow({col:row[col] for col in writer.fieldnames});match+=1
 assert seen==prev['common_valid_nonpalindromic'] and match==len(ids)
 results[k]={'source_sha256':sha(source),'source_rows':seen,'matched_rows':match,'output_path':str(out/f'{k}.sumstats.gz'),'output_sha256':sha(out/f'{k}.sumstats.gz')}
receipt={'selected_markers':len(ids),'ukb_hfrs_sha256':sha(raw),'sources':results,'ordered_snp_support_identical':True}
(p/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
