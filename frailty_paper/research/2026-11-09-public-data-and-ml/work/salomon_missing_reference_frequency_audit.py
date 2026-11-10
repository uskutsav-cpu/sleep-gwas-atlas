import csv,gzip,json,hashlib
from pathlib import Path
src=Path('work/salomon_neuronal_variant_readout.tsv.gz'); ld=Path('outputs/stage6_salomon_significant_region_ld.json');data=json.loads(ld.read_text());wanted={r['reporter_variant'] for r in data['results'] if r['reference_status']=='no_unique_biallelic_match'};found={}
with gzip.open(src,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  v=r['chr_pos_ref_alt']
  if v in wanted:
   o={k:r[k] for k in ['gnomad_AC','gnomad_AF','gnomad_AF_nfe','af_category','NGN2_variant_is_significant']}
   if v in found and found[v]!=o:raise ValueError('inconsistent source duplicate: '+v)
   found[v]=o
assert set(found)==wanted
out={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'ld_receipt_sha256':hashlib.sha256(ld.read_bytes()).hexdigest(),'count':len(wanted),'source_annotations':found,'limitations':'Annotations copied from processed source table and not separately verified in gnomAD. Absence of exact match in 1000G does not prove rarity in the target donors or establish allele error; no claim of disease effect or LD.'}
Path('outputs/stage6_salomon_missing_reference_frequency_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(found,indent=2))
