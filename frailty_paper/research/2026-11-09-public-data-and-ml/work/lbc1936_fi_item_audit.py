#!/usr/bin/env python3
import csv,hashlib,json,re
from pathlib import Path
from xml.etree import ElementTree as E
p=Path('work/lbc1936_inflammation_fi.xml');root=E.parse(p).getroot()
tables=[t for t in root.findall('.//table-wrap') if 'Items constituting the Frailty Index' in ''.join(t.itertext())]
assert len(tables)==1
rows=[]
for tr in tables[0].findall('.//tr'):
 cells=[' '.join(''.join(c.itertext()).split()) for c in tr if c.tag in ('th','td')]
 if len(cells)!=3 or cells[0]=='Items':continue
 item,coding,basis=cells
 rows.append({'item':item,'published_coding':coding,'basis':basis,'explicit_sleep_fatigue_term':bool(re.search(r'\b(sleep|insomnia|fatigue|tired|exhaustion|lethargy)\b',item,re.I)),'shared_physical_component_label':bool(re.search(r'grip|walk|chair|disability',item,re.I))})
assert len(rows)==30
with Path('outputs/lbc1936_published_fi_items.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
res={'source':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7456784/fullTextXML','xml_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'table':'A1','item_count':len(rows),'explicit_sleep_fatigue_item_labels':sum(r['explicit_sleep_fatigue_term'] for r in rows),'keyword_physical_function_labels':sum(r['shared_physical_component_label'] for r in rows),'note':'Keyword audit of top-level labels; underlying scales may contain related symptoms. Preserves source coding without resolving ambiguities.'}
Path('outputs/lbc1936_published_fi_items_receipt.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
