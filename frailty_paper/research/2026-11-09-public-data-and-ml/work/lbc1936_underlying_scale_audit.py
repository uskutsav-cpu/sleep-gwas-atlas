#!/usr/bin/env python3
import csv,hashlib,json,re,sys
sys.path.insert(0,'work/lbc_xls_reader')
import xlrd
p='work/lbc1936_public_dictionary.xls';s=xlrd.open_workbook(p).sheet_by_index(0)
rows=[]
for i in range(s.nrows):
 c=[str(x.value).strip() for x in s.row(i)]
 if re.fullmatch(r'HADS[AD][1-7]_w3',c[0]) or c[0] in {'hadsa_w3','hadsd_w3','hadstot_w3','adl_w3','chairst_w3','chairstone_w3'}:
  rows.append({'row':i+1,'variable':c[0],'label':c[1]})
items=[x for x in rows if re.fullmatch(r'HADS[AD][1-7]_w3',x['variable'])];assert len(items)==14
pattern=re.compile(r'\b(sleep|insomnia|fatigue|tired|exhaustion|lethargy)\b',re.I)
with open('outputs/lbc1936_fi_underlying_scale_fields.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['row','variable','label']);w.writeheader();w.writerows(rows)
summary={'source':'LBC1936 official Waves 1–4 public data dictionary','sha256':hashlib.sha256(open(p,'rb').read()).hexdigest(),'wave3_hads_item_count':len(items),'hads_items_explicit_keyword_matches':[x['variable'] for x in items if pattern.search(x['label'])],'fatigue_adjacent_item':'HADSD4_w3','townsend_individual_items_identified':False,'warning':'Labels alone do not determine complete constructs, coded responses, actual score composition, or the relationship between sleep and symptoms.'}
open('outputs/lbc1936_fi_underlying_scale_fields.json','w').write(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
