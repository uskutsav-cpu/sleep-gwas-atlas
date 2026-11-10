import csv,json,hashlib,math
from pathlib import Path
import pysam
root=Path('outputs');source=root/'stage6_aou_eur_candidate_association.csv'
rows=list(csv.DictReader(source.open()));candidates={}
for row in rows:
 if row['phenotype']!='NS_333.2':continue
 assert row['status']=='matched' and row['ancestry']=='eur' and row['sequencing_type']=='genome'
 rs=row['candidate'];assert rs not in candidates
 candidates[rs]={'chrom':row['chrom'],'pos':int(row['position']),'pair':{row['ref'],row['alt']}}
assert len(candidates)==33
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435387/harmonised/GCST90435387.h.tsv.gz'
# SSF header verified via bounded bytes from the same file before this lookup.
fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','Unnamed: 8','hm_coordinate_conversion','hm_code','variant_id']
t=pysam.TabixFile(url);out=[]
for rs,x in candidates.items():
 hits=[]
 try:raw=list(t.fetch(x['chrom'],x['pos']-1,x['pos']))
 except ValueError as e:raise RuntimeError('Missing target contig '+x['chrom']) from e
 for line in raw:
  vals=line.rstrip('\r\n').split('\t');assert len(vals)==len(fields);r=dict(zip(fields,vals))
  if int(r['base_pair_location'])==x['pos'] and {r['effect_allele'],r['other_allele']}==x['pair']:
   hits.append(r)
 if len(hits)>1:raise RuntimeError('Ambiguous exact pair '+rs)
 o={'variant':rs,'chrom_GRCh38':x['chrom'],'position_GRCh38':x['pos'],'expected_alleles':'/'.join(sorted(x['pair'])),'status':'no exact pair result'}
 if hits:
  r=hits[0];o.update(status='exact pair found',source_rsid=r['rsid'],effect_allele=r['effect_allele'],other_allele=r['other_allele'],beta=r['beta'],se=r['standard_error'],eaf=r['effect_allele_frequency'],p=r['p_value'],source_variant_id=r['variant_id'],hm_code=r['hm_code'])
  assert all(math.isfinite(float(o[k])) for k in ('beta','se','eaf','p')) and 0<float(o['p'])<=1
 out.append(o)
rank=sorted(range(33),key=lambda i:float(out[i].get('p',1)));prev=1.
for j in range(32,-1,-1):
 i=rank[j];prev=min(prev,float(out[i].get('p',1))*33/(j+1));out[i]['BH_q_33']=prev
fields=['variant','chrom_GRCh38','position_GRCh38','expected_alleles','status','source_rsid','effect_allele','other_allele','beta','se','eaf','p','BH_q_33','source_variant_id','hm_code']
p=root/'stage6_rls2024_candidate_lookup.csv'
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
receipt={'source':url,'catalog_accession':'GCST90435387','source_published_md5':'a220a11e2d14730988e561361d7cde94','source_whole_file_checksum_verified':False,'method':'Remote indexed range query; exact GRCh38 coordinate and unordered allele pair; output rsID supplied separately by source; missing assigned p=1 for exploratory BH over 33','candidate_input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'matches':sum(x['status']=='exact pair found' for x in out),'nominal_p_lt_0_05':sum(float(x.get('p',1))<.05 for x in out),'BH_q_lt_0_05':sum(x['BH_q_33']<.05 for x in out)}
(root/'stage6_rls2024_candidate_lookup_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
