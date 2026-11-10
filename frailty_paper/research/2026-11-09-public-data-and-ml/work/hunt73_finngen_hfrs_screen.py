#!/usr/bin/env python3
"""Descriptive rsID/allele screen; no directional or causal inference."""
import csv, hashlib, json, collections
from pathlib import Path
import openpyxl
root=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw')
s=openpyxl.load_workbook('work/hunt_circadian_supp.xlsx',read_only=True,data_only=True).active
hunt=collections.defaultdict(list)
for row in s.iter_rows(min_row=2,values_only=True):
 if isinstance(row[0],str) and row[0].startswith('rs'):hunt[row[0]].append(row)
assert len(hunt)==73 and all(len(x)==7 for x in hunt.values())
found=collections.defaultdict(list);h=hashlib.md5();n=0
with (root/'GCST90707291.tsv').open('rb') as f:
 header=None
 for line in f:
  h.update(line)
  if header is None:header=line.decode().strip().split('\t');continue
  n+=1
  vals=line.decode().strip().split('\t');rs=vals[-1]
  if rs in hunt: found[rs].append(dict(zip(header,vals)))
assert h.hexdigest()=='38d9434f50a4493def9d2c9d1cc1b4b3'
results=[]
for rs in sorted(hunt):
 a=hunt[rs][0]; source=found[rs]; pairs=[r for r in source if {r['effect_allele'].upper(),r['other_allele'].upper()}=={str(a[3]).upper(),str(a[4]).upper()}]
 for r in pairs or [{}]:
  results.append(dict(rsid=rs,hunt_gene=a[2],hunt_allele_pair=f'{a[3]}/{a[4]}',source_rows=len(source),pair_rows=len(pairs),hfrs_chr=r.get('chromosome',''),hfrs_pos_GRCh38=r.get('base_pair_location',''),hfrs_effect_allele=r.get('effect_allele',''),hfrs_other_allele=r.get('other_allele',''),hfrs_beta=r.get('beta',''),hfrs_se=r.get('standard_error',''),hfrs_p=r.get('p_value','')))
with open('outputs/stage6_hunt73_hfrs_rs_allele_screen.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
print(json.dumps({'hfrs_rows_scanned':n,'md5':h.hexdigest(),'unique_hunt_rsids':len(hunt),'rsids_seen':len(found),'with_pair_match':sum(bool(r['pair_rows']) for r in results),'pair_matched_p_lt_0_05':sum(bool(r['hfrs_p']) and float(r['hfrs_p'])<.05 for r in results),'multi_source_ids':sum(len(v)>1 for v in found.values())},indent=2))
