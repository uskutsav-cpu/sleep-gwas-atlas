import csv,hashlib,json,re
from pathlib import Path
import openpyxl
from pyliftover import LiftOver
out=Path(__file__).resolve().parent
rank=out.parents[1]/'prediction_rank_comparison/variant_ranks.tsv'
chain=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz')
book=Path('/Volumes/Extreme SSD/brain6-neuropsychiatric-mpra/data5.xlsx')
lo=LiftOver(str(chain)); data=list(csv.DictReader(rank.open(),delimiter='\t'))
w=openpyxl.load_workbook(book,read_only=True,data_only=True); raw={};summary=set()
for r in w['raw_counts'].iter_rows(values_only=True):
 if len(r)>1 and r[1]:
  m=re.match(r'^(chr\d+)\.(\d+)\.(Ref|all_alt)\.([ACGT])\.',str(r[1]))
  if m: raw.setdefault((m[1],int(m[2])),set()).add((m[3],m[4]))
for s in w:
 if s.title in ('raw_counts','correlation_samples','P-NPC'):continue
 for r in s.iter_rows(values_only=True):
  for x in r[:2]:
   m=re.fullmatch(r'(chr\d+)_(\d+)',str(x))
   if m:summary.add((m[1],int(m[2])))
rows=[]
for d in data:
 m=re.fullmatch(r'(chr\d+):(\d+):([ACGT])>([ACGT])',d['grch38'])
 if not m:
  rows.append({'rsid':d['snp'],'region':d['region'],'grch38':d['grch38'],'status':'invalid_snv'});continue
 chrom,pos,ref,alt=m.groups();lift=lo.convert_coordinate(chrom,int(pos)-1)
 if len(lift)!=1 or lift[0][0]!=chrom or lift[0][2]!='+':
  rows.append({'rsid':d['snp'],'region':d['region'],'grch38':d['grch38'],'status':'not_unique_forward_liftover'});continue
 pos37=lift[0][1]+1; key=(chrom,pos37); alleles=raw.get(key,set())
 status='matched_ref_alt' if {('Ref',ref),('all_alt',alt)} <= alleles else ('coordinate_only_allele_mismatch' if alleles else 'no_raw_construct_match')
 rows.append({'rsid':d['snp'],'region':d['region'],'grch38':d['grch38'],'grch37':f'{chrom}:{pos37}','status':status,'summary_coordinate_present':key in summary,'construct_alleles':';'.join(':'.join(a) for a in sorted(alleles))})
fields=['rsid','region','grch38','grch37','status','summary_coordinate_present','construct_alleles']
with (out/'ranked_screen.tsv').open('w',newline='') as f:
 dw=csv.DictWriter(f,fieldnames=fields,delimiter='\t');dw.writeheader();dw.writerows(rows)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
from collections import Counter
receipt={'rank_sha256':sha(rank),'chain_sha256':sha(chain),'workbook_sha256':sha(book),'source_raw_coordinates_with_named_alleles':len(raw),'status_counts':dict(Counter(r['status'] for r in rows)),'summary_sheets_excluded':['P-NPC','correlation_samples'],'coordinate_system':'input GRCh38 one-based; chain conversion zero-based; output GRCh37 one-based','source_note':'MPRA genomic oligo design described as hg37 in paper; reference and alternative identities matched by raw construct label','output_sha256':sha(out/'ranked_screen.tsv')}
(out/'ranked_screen_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
