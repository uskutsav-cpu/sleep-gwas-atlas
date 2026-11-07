#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parents[2]
source=Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2');loc='placo_insomnia__adhd_chr3_117142005_118142005'
prov=json.loads((root.parent/'brain6-alphagenome-pilot/results/provenance.json').read_text())
ref=json.loads((p.parent/'finngen_supplement/reference_ld/receipt.json').read_text())
def sha(x):
 h=hashlib.sha256()
 with x.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
variants=ref['variants'];assert len(variants)==6
pipfile=source/'variant_pip.tsv';assert sha(pipfile)==prov['full_variant_pip_sha256']
pip={}
with pipfile.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['candidate_locus_id']==loc and r['snp'] in {v['rsid'] for v in variants.values()}:
   k=(r['snp'],r['trait_id']);assert k not in pip;pip[k]=float(r['conditional_pip'])
results=[];inputs={}
for trait in ('insomnia','adhd'):
 filename=loc+'__'+trait+'.tsv';path=source/'inputs'/filename;assert sha(path)==prov['input_sha256'][filename]
 inputs[trait]={'path':str(path),'sha256':sha(path)}
 index={}
 with path.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['SNP'] in {v['rsid'] for v in variants.values()}:index.setdefault(r['SNP'],[]).append(r)
 for label,v in variants.items():
  matches=[r for r in index.get(v['rsid'],[]) if r['CHR']=='3' and int(r['POS'])==v['pos37'] and {r['A1'],r['A2']}=={v['ref'],v['alt']}]
  assert len(matches)==1 and (v['rsid'],trait) in pip
  x=matches[0];beta,se=float(x['BETA']),float(x['SE']);assert se>0
  results.append({'variant_role':label,'snp':v['rsid'],'trait':trait,'pos37':v['pos37'],
   'a1':x['A1'],'a2':x['A2'],'beta_a1':beta,'se':se,'wald_p_approx':math.erfc(abs(beta/se)/math.sqrt(2)),
   'conditional_pip':pip[(v['rsid'],trait)]})
assert len(results)==12
with (p/'comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=results[0].keys(),delimiter='\t');w.writeheader();w.writerows(results)
receipt={'input_files':inputs,'pip_path':str(pipfile),'pip_sha256':sha(pipfile),'ld_receipt_sha256':sha(p.parent/'finngen_supplement/reference_ld/receipt.json'),
 'comparison_sha256':sha(p/'comparison.tsv'),'matched_rows':len(results),'note':'Selected markers; PIP conditional on stored single-signal model; Wald p approximated from beta/SE.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
