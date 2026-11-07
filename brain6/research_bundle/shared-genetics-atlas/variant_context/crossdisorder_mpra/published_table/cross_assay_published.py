#!/usr/bin/env python3
import csv,hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
old=read(p.parent/'cross_assay_candidate_comparison.tsv');lee={x['rsid']:x for x in read(p/'candidate_comparison.tsv')}
assert len(old)==4 and len({x['rsid'] for x in old})==4
out=[]
for x in old:
 r=lee[x['rsid']];m=re.fullmatch(r'(chr\d+):(\d+)_([ACGT]+):([ACGT]+)',r['source_key']);assert m
 assert x['variant37']==':'.join(m.groups()) and r['allelic_status']=='unique_match' and r['element_status']=='unique_match'
 out.append({'rsid':x['rsid'],'variant37':x['variant37'],'lee_allelic_fdr':r['allelic_fdr'],'lee_variant_class':r['allelic_class'],'lee_element_fdr':r['element_fdr'],'lee_element_class':r['element_class'],'siraj_library':x['siraj_library'],'siraj_allelic_fdr':x['siraj_source_fdr'],'siraj_a_activity_p':x['siraj_A_activity_bonferroni_p'],'siraj_b_activity_p':x['siraj_B_activity_bonferroni_p']})
with (p/'cross_assay_published.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
both=[x['rsid'] for x in out if float(x['lee_allelic_fdr'])<.05 and float(x['siraj_allelic_fdr'])<.05]
receipt={'candidate_count':len(out),'both_reported_allelic_fdr_lt_005':both,'rows_sha256':hashlib.sha256((p/'cross_assay_published.tsv').read_bytes()).hexdigest(),'limitations':'Post-selected, different assay contexts and testing families; allele-effect sign and authors activity criteria are not reconciled; not independent genetic replication.'}
(p/'cross_assay_published_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
