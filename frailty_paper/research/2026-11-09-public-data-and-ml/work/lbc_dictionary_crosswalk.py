#!/usr/bin/env python3
import csv,hashlib,json,sys
sys.path.insert(0,'work/lbc_xls_reader')
import xlrd
p='work/lbc1936_public_dictionary.xls'; w=xlrd.open_workbook(p); s=w.sheet_by_index(0)
targets={'sqquality_w3','sqquicklywknights_w3','sqquicklywkends_w3','sqawakening_w3','sqlengthwknights_w3','sixmwk_w3','griprh_1_w3','griprh_2_w3','griprh_3_w3','griplh_1_w3','griplh_2_w3','griplh_3_w3','phyactiv_w3','sqawakening_excToilet_w4','sqawakening_specToilet_w4','phyactiv_w4'}
rows=[]
for i in range(s.nrows):
 r=[str(c.value).strip() for c in s.row(i)]
 if r[0] in targets: rows.append({'dictionary_row':i+1,'variable':r[0],'description':r[1],'type':r[2]})
assert len(rows)==len(targets) and len({r['variable'] for r in rows})==len(targets)
with open('outputs/lbc1936_dictionary_wave_crosswalk.csv','w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
print(json.dumps({'rows':len(rows),'dictionary_sha256':hashlib.sha256(open(p,'rb').read()).hexdigest(),'sheets':w.sheet_names()},indent=2))
