"""Bounded read-only XLSX extraction with bundled Python; no workbook edits.

Streaming XML avoids importing millions of spreadsheet cells into another
in-memory workbook while the measured-effect model is running.
"""
import csv,json,hashlib,posixpath,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
BASE=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v7')
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def inspect(path,max_rows=8,print_rows=5):
    result={'file':str(path),'sha256':digest(path),'sheets':[]}
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None;strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for event,node in ET.iterparse(z.open('xl/sharedStrings.xml'),events=('end',)):
                if node.tag.endswith('}si'):strings.append(''.join(node.itertext()));node.clear()
        relationships={n.attrib['Id']:n.attrib['Target'] for n in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};book=ET.fromstring(z.read('xl/workbook.xml'))
        for sheet in book.find('m:sheets',NS):
            rid=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'];target=relationships[rid];target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target);rows=[];dimension=None
            for event,node in ET.iterparse(z.open(target),events=('end',)):
                if node.tag.endswith('}dimension'):dimension=node.attrib.get('ref')
                if not node.tag.endswith('}row'):continue
                cells=[]
                for cell in node:
                    if not cell.tag.endswith('}c'):continue
                    value=cell.find('m:v',NS);text=value.text if value is not None else '';kind=cell.attrib.get('t','n')
                    if kind=='s':text=strings[int(text)] if text else ''
                    elif kind=='inlineStr':text=''.join(cell.find('m:is',NS).itertext())
                    formula=cell.find('m:f',NS);cells.append({'cell':cell.attrib.get('r'),'value':text,'type':kind,**({'formula':formula.text} if formula is not None else {})})
                rows.append({'row':node.attrib.get('r'),'cells':cells});node.clear()
                if len(rows)>=max_rows:break
            item={'name':sheet.attrib['name'],'dimension':dimension,'first_rows':rows};result['sheets'].append(item);print(path.name,'SHEET',item['name'],'dimension',dimension,flush=True)
            for row in rows[:print_rows]:print(' ROW',row['row'],[(c['cell'],str(c['value'])[:170]) for c in row['cells'][:30]],flush=True)
    return result
def main():
    names=['S3X5WEP33_NIHMS787218-supplement-7_TableS1.xlsx','S3X5WEP33_NIHMS787218-supplement-9_TableS3.xlsx','S12C1A15T_PMC11839127_media-3_ENIGMA_MPRA_results.xlsx'];reports=[inspect(BASE/'sources'/n) for n in names];p=BASE/'outputs/additional_source_workbook_structure_audit.json';p.write_text(json.dumps({'files':reports,'method':'Read-only streaming XLSX XML extraction with source ZIP CRC checks; original spreadsheets untouched.','scientific_status':'Source structure only, not an accepted validation dataset or modeled performance result.'},indent=2)+'\n')
if __name__=='__main__':main()
