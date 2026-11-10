import csv,json,hashlib,re
from pathlib import Path
p=Path('outputs');maps={}
for fn,cols in [('integrated_15_candidate_review.csv',('GRCh37_chr','GRCh38_pos')),('expanded_22_candidate_risk_gene_table.csv',('GRCh38_chr','GRCh38_pos'))]:
 for r in csv.DictReader((p/fn).open()):
  key=(r[cols[0]],r[cols[1]],frozenset((r['counted_allele'],r['other_allele'])));variant=r['variant']
  if variant in maps:assert maps[variant]==key
  maps[variant]=key
assert len(maps)==33 and len(set(maps.values()))==33
coverage=p/'stage5_siraj2026_allcell_coverage.csv'; lookup={}
for r in csv.DictReader(coverage.open()):
 m=re.fullmatch(r'chr(\d+):(\d+):([ACGT]+)>([ACGT]+)',r['variant_GRCh38'])
 if not m:continue
 key=(m[1],m[2],frozenset((m[3],m[4])),r['cell_type']);assert key not in lookup;lookup[key]=r
out=[]
for variant,(ch,bp,pair) in maps.items():
 for cell in ['SKNSH','K562','HEPG2','HCT116','A549','GM12878']:
  r=lookup.get((ch,bp,pair,cell)); out.append({'variant':variant,'GRCh38_chr':ch,'GRCh38_pos':bp,'alleles_unordered':'/'.join(sorted(pair)),'cell_type':cell,'coverage_status':r['grch37_mapping_status'] if r else 'not present in queue coverage', 'selected_record_groups':r['selected_record_groups'] if r else '', 'author_rule_classification':r['selected_emvar'] if r else 'not assessed','source_variant':r['variant_GRCh38'] if r else '', 'caveat':'Non-emVar is not a validated negative; no record is missing not zero; selected reporter assay is not clinical validation.'})
with (p/'stage6_33_candidate_reporter_coverage.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
from collections import Counter
res={'variant_count':len(maps),'variant_cell_rows':len(out),'counts':dict(Counter(x['author_rule_classification'] for x in out)),'inputs':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ['integrated_15_candidate_review.csv','expanded_22_candidate_risk_gene_table.csv',coverage.name]},'matching':'GRCh38 position and unordered exact allele pair; no proxies; study selection as documented in source coverage'}
(p/'stage6_33_candidate_reporter_coverage.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
