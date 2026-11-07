#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent
ld=json.loads((p.parent/'finngen_supplement/reference_ld/receipt.json').read_text())['variants']
local=p.parent/'fine_mapping_comparison/comparison.tsv'
with local.open() as f:brain={r['snp']:r for r in csv.DictReader(f,delimiter='\t') if r['trait']=='adhd'}
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90568001-GCST90569000/GCST90568440/harmonised/GCST90568440.h.tsv.gz'
cols=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','hm_coordinate_conversion','hm_code','variant_id']
tab=pysam.TabixFile(url);out=[];queries=[]
for label,v in ld.items():
 original=brain[v['rsid']];lines=list(tab.fetch('3',v['pos38']-1,v['pos38']));parsed=[]
 for line in lines:
  z=line.split('\t');assert len(z)==len(cols);parsed.append(dict(zip(cols,z)))
 hits=[r for r in parsed if r['rsid']==v['rsid'] and int(r['base_pair_location'])==v['pos38'] and {r['effect_allele'],r['other_allele']}=={v['ref'],v['alt']}]
 row={'role':label,'snp':v['rsid'],'pos38':v['pos38'],'brain_a1':original['a1'],'brain_a2':original['a2'],
  'brain_adhd_beta_a1':original['beta_a1'],'status':'MISSING_OR_AMBIGUOUS','symptom_beta_a1':'','symptom_se':'','symptom_p':'','direction_agreement':''}
 if len(hits)==1:
  h=hits[0];b,se,pval=map(float,(h['beta'],h['standard_error'],h['p_value']))
  assert all(map(math.isfinite,(b,se,pval))) and se>0 and 0<=pval<=1
  b*=1 if h['effect_allele']==original['a1'] else -1
  row.update(status='MATCHED',symptom_beta_a1=b,symptom_se=se,symptom_p=pval,
   direction_agreement='YES' if b*float(original['beta_a1'])>0 else 'NO_OR_ZERO')
 out.append(row);queries.append({'snp':v['rsid'],'raw_lines':lines,'unique_matches':len(hits)})
with (p/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t');w.writeheader();w.writerows(out)
receipt={'url':url,'query_records':queries,'table_sha256':hashlib.sha256((p/'results.tsv').read_bytes()).hexdigest(),
 'brain_comparison_sha256':hashlib.sha256(local.read_bytes()).hexdigest(),
 'source_note':'Remote tabix and index not whole-file hashed; source cohort independence is author-reported.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(out,indent=2))
