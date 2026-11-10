#!/usr/bin/env python3
"""Match exact queue alleles to public Ho et al. 2025 supplemental tables."""
import argparse,csv,openpyxl,hashlib,json
from collections import defaultdict,Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--workbook',required=True);a=p.parse_args();wb=Path(a.workbook)
with open('outputs/stage5_score_queue.csv',newline='') as f:queue=list(csv.DictReader(f))
qidx={(x['chr'],int(x['position_GRCh38']),x['reference_allele'],x['alternate_allele']):x for x in queue};by_pos=defaultdict(list)
for k in qidx:by_pos[k[:2]].append(k)
w=openpyxl.load_workbook(wb,read_only=True,data_only=True);out=[];positions=defaultdict(int);summary={}
for sheet in ['1. Tcell MPRA results','9. jurkat MPRA results']:
 s=w[sheet];it=s.iter_rows(values_only=True);header=next(it);assert len(header)==len(set(header));n=0
 for row in it:
  d=dict(zip(header,row));n+=1
  try:k=(str(d['chr']),int(d['pos']),str(d['ref_allele']).upper(),str(d['alt_allele']).upper())
  except (KeyError,TypeError,ValueError):continue
  if k[:2] in by_pos:positions[(sheet,k[:2])]+=1
  q=qidx.get(k)
  if not q:continue
  out.append({'variant_key':q['variant_key'],'source_sheet':sheet,'source_row':n+1,'rsid':d.get('rsid'),'SNP19':d.get('SNP19'),'SNP38':d.get('SNP38'),'ref_allele':d.get('ref_allele'),'alt_allele':d.get('alt_allele'),'mpra_sig':d.get('mpra_sig'),'strand':d.get('strand'),'allele':d.get('allele'),'haplotype':d.get('haplotype'),'A_log2FC':d.get('A_log2FC'),'B_log2FC':d.get('B_log2FC'),'Log2Skew':d.get('Log2Skew'),'Skew_logFDR':d.get('Skew_logFDR'),'Skew_SE':d.get('Skew_SE'),'Skew_logP':d.get('Skew_logP'),'interpretation':'Exact GRCh38 allele labels; experimentally measured in reporter; verify direction and independence before sign comparison'})
 summary[sheet]={'source_rows':n,'exact_matched_records':sum(x['source_sheet']==sheet for x in out),'exact_variants':len({x['variant_key'] for x in out if x['source_sheet']==sheet})}
fields=['variant_key','source_sheet','source_row','rsid','SNP19','SNP38','ref_allele','alt_allele','mpra_sig','strand','allele','haplotype','A_log2FC','B_log2FC','Log2Skew','Skew_logFDR','Skew_SE','Skew_logP','interpretation']
with open('outputs/stage5_ho2025_mpra_exact_alleles.csv','w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(out)
print(json.dumps(summary,indent=2));print('categories',dict(Counter((r['source_sheet'],r['mpra_sig']) for r in out)))
