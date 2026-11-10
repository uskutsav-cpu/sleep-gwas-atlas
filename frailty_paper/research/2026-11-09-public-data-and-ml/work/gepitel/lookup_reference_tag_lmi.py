#!/usr/bin/env python3
"""Look up the previous selected tags in the released LMI workbook."""
import csv,hashlib,json
from pathlib import Path
import openpyxl
source=Path('work/gepitel/GWAS_LMI.xlsx'); tags=Path('outputs/stage6_gepitel_reference_tag_call_summary.csv')
h=hashlib.md5()
with source.open('rb') as f:
 for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
assert h.hexdigest()=='c58c6d9b4427fc9b16c223fe4bc1faae'
rows=list(csv.DictReader(tags.open())); target={r['reference_tag'] for r in rows}; found={k:[] for k in target}
s=openpyxl.load_workbook(source,read_only=True,data_only=True)['dataset'];n=0; p_valid=0; p_ge05=0; max_p=0.0; p_invalid=0
for row in s.iter_rows(min_row=2,values_only=True):
 n+=1
 try:
  pval=float(row[6]); assert 0<=pval<=1; p_valid+=1; p_ge05+=pval>=0.05; max_p=max(max_p,pval)
 except (TypeError,ValueError,AssertionError): p_invalid+=1
 if row[1] in found:
  found[row[1]].append(dict(zip(['source_chr','source_rsid','source_position','source_alleleA','source_alleleB','source_info','source_p','source_beta','source_se'],row)))
out=[]
for r in rows:
 f=found[r['reference_tag']]
 out.append({'candidate':r['candidate'],'score_lists':r['score_lists'],'reference_tag':r['reference_tag'],'reference_r2':r['reference_r2'],'source_id_matches':len(f),**({k:v for k,v in f[0].items()} if len(f)==1 else {})})
fields=['candidate','score_lists','reference_tag','reference_r2','source_id_matches','source_chr','source_rsid','source_position','source_alleleA','source_alleleB','source_info','source_p','source_beta','source_se']
with open('outputs/stage6_gepitel_tag_lmi_summary.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
receipt={'workbook_md5':h.hexdigest(),'sheet':'dataset','rows_scanned':n,'p_valid':p_valid,'p_invalid':p_invalid,'p_ge_0_05':p_ge05,'max_p':max_p,'selected_tag_ids':len(target),'tag_ids_found':sum(bool(v) for v in found.values()),'tag_ids_duplicate':sum(len(v)>1 for v in found.values()),'tag_list_sha256':hashlib.sha256(tags.read_bytes()).hexdigest(),'note':'rsID identifier match only; source build and effect allele unresolved; reference LD does not establish candidate or score effect.'}
Path('outputs/stage6_gepitel_tag_lmi_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps(receipt,indent=2))
