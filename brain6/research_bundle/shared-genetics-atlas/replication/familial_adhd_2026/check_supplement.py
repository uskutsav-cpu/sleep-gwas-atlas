from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import hashlib,json
p=Path(__file__).resolve().parent;f=p/'source_supplement.docx';ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with ZipFile(f) as z:
 bad=z.testzip();assert bad is None;root=ET.fromstring(z.read('word/document.xml'))
text=' '.join(x.text or '' for x in root.findall('.//w:t',ns));queries=['rs62264764','rs62264767','rs2431108','rs77960','rs1476535','rs2582894','rs7160630','rs6452785','LSAMP','FOXP2','MEF2C'];counts={x:text.lower().count(x.lower()) for x in queries}
tables=[]
for i,t in enumerate(root.findall('.//w:tbl',ns),1):
 rows=[]
 for row in t.findall('./w:tr',ns):rows.append([' '.join(a.text or '' for a in c.findall('.//w:t',ns)) for c in row.findall('./w:tc',ns)])
 hits=[]
 for j,row in enumerate(rows,1):
  if any('foxp2' in x.lower() or 'mef2c' in x.lower() for x in row):hits.append({'row':j,'first_cell':row[0] if row else ''})
 tables.append({'table':i,'rows':len(rows),'gene_name_hits':hits})
out={'source':'https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs00439-026-02840-7/MediaObjects/439_2026_2840_MOESM1_ESM.docx','sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size,'zip_integrity':'passed','text_substring_counts':counts,'tables':tables,'scope':'Word document text only; images, external datasets, exact per-variant p-values and cohort overlap not established.'};(p/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(counts,indent=2))
