#!/usr/bin/env python3
import collections,csv,gzip,hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent;source=p/'cre_sequences.tsv.gz';rank=p.parents[1]/'prediction_rank_comparison/variant_ranks.tsv';region=p.parent/'rs62264764_dense_qtl/multisignal_pilot/variants.tsv'
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
with rank.open() as f:candidates=list(csv.DictReader(f,delimiter='\t'))
ids={x['snp'] for x in candidates};ids.add('rs62264767')
with region.open() as f:regional=list(csv.DictReader(f,delimiter='\t'))
regional_ids={x['rsid'] for x in regional};assert len(regional)==len(regional_ids)==1548
counts=collections.Counter();hits=collections.defaultdict(set);regional_hits=collections.defaultdict(set);unexpected=[];records=0
with gzip.open(source,'rt') as f:
 reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==['names','seq']
 for row in reader:
  records+=1;name=row['names'];m=re.fullmatch(r'(rs\d+)_scz_(ref|alt)',name)
  if not m:unexpected.append(name);continue
  rs,allele=m.groups();counts[allele]+=1
  if rs in ids:hits[rs].add(allele)
  if rs in regional_ids:regional_hits[rs].add(allele)
report={'source_url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211045/suppl/GSE211045_cre_sequences.tsv.gz','compressed_bytes':source.stat().st_size,'source_sha256':sha(source),'candidate_file_sha256':sha(rank),'candidate_rows':len(candidates),'unique_candidate_rsids_plus_neighbor':len(ids),'regional_file_sha256':sha(region),'regional_unique_rsids':len(regional_ids),'construct_rows':records,'allele_label_counts':dict(counts),'unexpected_name_count':len(unexpected),'unexpected_name_examples':unexpected[:10],'candidate_matches':{x:sorted(y) for x,y in sorted(hits.items())},'regional_matches':{x:sorted(y) for x,y in sorted(regional_hits.items())},'limitations':'Name-only intersection; unknown build/allele orientation, testing/QC status, and biological effects. Complete gzip read checks stream CRC; no external publisher hash was available. Regional extension was post hoc.'}
(p/'lookup_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
