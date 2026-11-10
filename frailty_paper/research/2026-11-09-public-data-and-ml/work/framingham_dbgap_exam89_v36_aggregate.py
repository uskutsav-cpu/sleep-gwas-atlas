#!/usr/bin/env python3
import csv, xml.etree.ElementTree as ET
sources=[('8','work/framingham_dbgap_exam8_v36_report.xml',{'H401','H471','H472','H480','H481','H482','H483','H484','H605','H617','H730','H731','H732','H745'}),('9/Omni 4','work/framingham_dbgap_exam9_v36_report.xml',{'J476','J628','J629','J630','J631','J632','J733','J746','J765','J766','J767','J768','J769','J770','J785','J790','J944'})]
rows=[]
for exam,path,names in sources:
 root=ET.parse(path).getroot(); seen=set()
 for v in root.findall('variable'):
  name=v.get('var_name','').upper()
  if name not in names or '.c' in v.get('id',''): continue
  stats=v.find('total/stats'); stat=stats.find('stat');seen.add(name)
  rows.append({'exam':exam,'dataset':root.get('dataset_id'),'variable':name,'accession':v.get('id'),'description':v.findtext('description'),'nonmissing':stat.get('n'),'missing':stat.get('nulls'),'category_counts':'; '.join(f"{x.get('code')}={x.get('count')}" for x in stats.findall('enum'))})
 assert seen==names,(path,names-seen)
with open('outputs/stage6_framingham_exam89_dbgap_v36_aggregate.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
