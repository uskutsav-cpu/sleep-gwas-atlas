#!/usr/bin/env python3
from pathlib import Path
import requests,hashlib,json,pysam,csv
O=Path('outputs');W=Path('work');i=requests.get('https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/7134291',timeout=30);i.raise_for_status();j=i.json();placements=[p for p in j['primary_snapshot_data']['placements_with_allele'] if p.get('is_ptlp') and p['seq_id']=='NC_000012.12'];assert len(placements)==1;a=placements[0]['alleles'];vals=[x['allele']['spdi'] for x in a];pos={x['position']+1 for x in vals};assert len(pos)==1;pos=pos.pop();alleles={x['inserted_sequence'] for x in vals}|{x['deleted_sequence'] for x in vals};print('NCBI GRCh38',pos,sorted(alleles),flush=True);assert {'A','G'}<=alleles
spec=[('GCST90475826','GCST90475001-GCST90476000','MVP insomnia'),('GCST90707291','GCST90707001-GCST90708000','FinnGen hospital frailty')];results=[];receipts=[]
for acc,group,trait in spec:
 base=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{group}/{acc}/harmonised/';name=f'{acc}.h.tsv.gz';r=requests.get(base+name+'.tbi',timeout=45);r.raise_for_status();idx=W/(name+'.tbi');idx.write_bytes(r.content);m=requests.get(base+'md5sum.txt',timeout=30);m.raise_for_status();md5=hashlib.md5(r.content).hexdigest();expected=[line.split()[0] for line in m.text.splitlines() if len(line.split())==2 and line.split()[-1]==name+'.tbi'];assert not expected or expected==[md5];f=pysam.TabixFile(base+name,index=str(idx));lines=list(f.fetch('12',pos-1,pos));cands=[]
 for line in lines:
  x=line.split('\t');
  if x[1]==str(pos) and {x[2],x[3]}=={'A','G'}:cands.append(x)
 row={'trait':trait,'accession':acc,'rsid':'rs7134291','chr38':12,'pos38':pos,'requested_counted_allele':'A','other_allele':'G','position_rows':len(lines),'exact_allele_rows':len(cands)}
 if len(cands)==1:
  x=cands[0];row.update(source_effect_allele=x[2],source_other_allele=x[3],source_stat=x[4],source_se=x[5],source_eaf=x[6],p_value=x[7],trailing_fields='|'.join(x[8:]));
  # Preserve raw effect; OR vs beta is determined only after original metadata/header audit.
 results.append(row);receipts.append({'accession':acc,'base':base,'index_md5':md5,'index_matches_published_md5': (expected==[md5]) if expected else None,'full_file_checksum_verified':False});print(acc,row,flush=True)
fields=['trait','accession','rsid','chr38','pos38','requested_counted_allele','other_allele','position_rows','exact_allele_rows','source_effect_allele','source_other_allele','source_stat','source_se','source_eaf','p_value','trailing_fields']
with (O/'stage6_hrs_elsa_lead_external_lookup.csv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(results)
(O/'stage6_hrs_elsa_lead_external_receipt.json').write_text(json.dumps({'variant_ncbi_url':i.url,'variant_ncbi_sha256':hashlib.sha256(i.content).hexdigest(),'source_receipts':receipts,'coordinate_system':'GRCh38 1-based','selected_locus':'Mekli et al. HRS/ELSA FI nominally replicated rs7134291','caution':'postselected sleep lookup; HRS/ELSA sample overlap with other sources not established; exact score not tested'},indent=2)+'\n')
