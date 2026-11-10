import csv,json,hashlib,math
from pathlib import Path
import pysam
O=Path('outputs');p=O/'stage6_rls2024_candidate_lookup.csv';targets={r['variant']:r for r in csv.DictReader(p.open()) if r['variant'] in ['rs113851554','rs78206187']};assert len(targets)==2
fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','hm_coordinate_conversion','hm_code','variant_id']
sets=[('GCST90399568','EU-RLS-GENE'),('GCST90399571','INTERVAL')];rows=[];sources={}
for acc,group in sets:
 url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90399001-GCST90400000/{acc}/harmonised/{acc}.h.tsv.gz';sources[acc]=url;t=pysam.TabixFile(url)
 for rs,x in targets.items():
  candidates=[]
  for ln in t.fetch(x['chrom_GRCh38'],int(x['position_GRCh38'])-1,int(x['position_GRCh38'])):
   v=ln.rstrip('\r\n').split('\t');assert len(v)==len(fields);v=dict(zip(fields,v))
   if {v['effect_allele'],v['other_allele']}==set(x['expected_alleles'].split('/')):candidates.append(v)
  if len(candidates)>1:raise RuntimeError('ambiguous '+rs)
  o={'study':acc,'cohort':group,'variant':rs,'status':'no exact pair'}
  if candidates:
   v=candidates[0]; assert v['rsid']==rs and v['effect_allele'] in (x['effect_allele'],x['other_allele']); beta=float(v['beta'])*(1 if v['effect_allele']==x['effect_allele'] else -1)
   o.update(status='exact pair',counted_allele=x['effect_allele'],study_effect_allele=v['effect_allele'],study_other_allele=v['other_allele'],counted_beta=beta,se=v['standard_error'],p=v['p_value'],effect_allele_frequency=v['effect_allele_frequency'],variant_id=v['variant_id'],hm_code=v['hm_code']); assert 0<float(o['p'])<=1 and math.isfinite(beta)
  rows.append(o)
fields=['study','cohort','variant','status','counted_allele','study_effect_allele','study_other_allele','counted_beta','se','p','effect_allele_frequency','variant_id','hm_code'];out=O/'stage6_rls_schormair_two_variants.csv'
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
receipt={'sources':sources,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'source_full_file_checksums_verified':False,'description':'Two candidate lookup in two CC0 source cohorts, position and allele matched and beta oriented to prototype allele; tabix remote ranges','caveats':'Sleep movement associations only; participant overlap uncertain; no individual score, no causal or predictive inference; selected targets.'};(O/'stage6_rls_schormair_two_variants_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(rows,indent=2))
