import xml.etree.ElementTree as ET,json,hashlib
from pathlib import Path
p=Path('work/adni_fi_paper.xml');root=ET.parse(p).getroot(); found=[]
for t in root.findall('.//table-wrap'):
 rows=[[' '.join(' '.join(x.itertext()).split()) for x in tr] for tr in t.findall('.//tr')]
 if any(r and r[0]=='Self-reported insomnia' for r in rows):found.append(rows)
if len(found)!=1:raise ValueError('Unexpected number of FI item tables')
rows=found[0];wanted=['NPI, sleep disturbance','Self-reported drowsiness','Self-reported low energy','Self-reported insomnia'];result={}
for label in wanted:
 matches=[r for r in rows if r and r[0]==label]
 if len(matches)!=1 or len(matches[0])!=4:raise ValueError('Unexpected row '+label)
 result[label]=dict(zip(['FI_93','FI_26_refined','FI_40'],[x=='X' for x in matches[0][1:]]))
if sum(len(r)==4 and r[2]=='X' for r in rows)!=26:raise ValueError('26-item index mismatch')
out={'source':'https://doi.org/10.1007/s11357-022-00669-2','source_xml_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'items':result,'limits':'Paper-defined indices only; no participant data or clinical score validation. NPI and NPI-Q sleep domains may overlap; exact phase-specific mappings require verification.'}
Path('outputs/stage6_adni_fi_sleep_item_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(result,indent=2))
