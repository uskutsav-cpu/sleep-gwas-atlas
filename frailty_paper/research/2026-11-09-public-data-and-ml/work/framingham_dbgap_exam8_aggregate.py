#!/usr/bin/env python3
"""Extract non-consent-stratified aggregate metadata; no individual data."""
import csv, xml.etree.ElementTree as ET
from pathlib import Path
p=Path('work/framingham_dbgap_exam8_v20_report.xml')
r=ET.parse(p).getroot()
names={'H401','H605','H617','H730','H731','H732','H745'}
rows=[]
for v in r.findall('variable'):
 name=v.get('var_name','').upper()
 if name not in names or '.c' in v.get('id',''): continue
 stats=v.find('total/stats'); stat=stats.find('stat')
 rows.append({'variable':name,'dbgap_variable':v.get('id'),'description':v.findtext('description'),'nonmissing':stat.get('n'),'missing':stat.get('nulls'),'category_counts':'; '.join(f"{x.get('code')}={x.get('count')}" for x in stats.findall('enum'))})
assert {x['variable'] for x in rows} == names
with open('outputs/stage6_framingham_exam8_dbgap_aggregate.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
