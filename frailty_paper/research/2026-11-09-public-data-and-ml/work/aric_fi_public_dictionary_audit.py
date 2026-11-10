import csv, json, hashlib, re, xml.etree.ElementTree as ET
from pathlib import Path
p=Path('work/aric_findexcomp51_public_dict.xml')
root=ET.parse(p).getroot()
rows=[]
for v in root.findall('variable'):
 name=v.findtext('name','');desc=v.findtext('description','')
 if re.fullmatch(r'FDI[0-9]+',name):
  rows.append({'variable':name,'accession':v.attrib['id'],'description':desc,'sleep_keyword':bool(re.search(r'sleep|insomnia|restless|circadian',desc,re.I)),'exhaustion_question':name in {'FDI10','FDI14'}})
assert len(rows)==59 and {x['variable'] for x in rows}=={f'FDI{i}' for i in range(1,60)}
with open('outputs/aric_fi_public_item_inventory.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({'items':len(rows),'sleep_keyword_matches':[x['variable'] for x in rows if x['sleep_keyword']],'exhaustion_questions':[x['variable'] for x in rows if x['exhaustion_question']],'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
